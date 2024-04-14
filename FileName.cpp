#define _CRT_SECURE_NO_WARNINGS
#include <algorithm>
#include <atomic>
#include <chrono>
#include <iostream>
#include <ostream>
#include <sstream>
#include <string>
#include <random>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <fstream>
#include <cstring>
#include <cstdint>
#include <cassert>

using namespace std::chrono_literals;


constexpr int64_t operator""_kb(unsigned long long value) {
    return value * 1024;
}

constexpr int64_t operator""_mb(unsigned long long value) {
    return value * (1024 * 1024);
}

constexpr int64_t operator""_gb(unsigned long long value) {
    return value * (1024 * 1024 * 1024);
}

/// Allocate aligned for @count objects of type T, does not perform initialization
/// @param count - the number of objects
/// @param unaligned [out] - stores the un-aligned pointer, used to call free
/// @return pointer to the memory or nullptr
template <typename T>
T *alignedAlloc(size_t count, void *&unaligned) {
    const size_t bytes = count * sizeof(T);
    unaligned = malloc(bytes + 63);
    if (!unaligned) {
        return nullptr;
    }
    T *const aligned = reinterpret_cast<T *>(uintptr_t(unaligned) + 63 & -64);
    return aligned;
}

/// Stack allocator with predefined max size
/// The total memory is 64 byte aligned, all but the first allocation are not guaranteed to be algigned
/// Can only free all the allocations at once
struct StackAllocator {
    StackAllocator(uint8_t *ptr, int64_t bytes)
        : totalBytes(bytes)
        , data(ptr) {
    }

    /// Allocate memory for @count T objects
    /// Does *NOT* call constructors
    /// @param count - the number of objects needed
    /// @return pointer to the allocated memory or nullptr
    template <typename T>
    T *alloc(int64_t count) {
        const int64_t size = count * sizeof(T);
        if (idx + size > totalBytes) {
            return nullptr;
        }
        uint8_t *start = data + idx;
        idx += size;
        return reinterpret_cast<T *>(start);
    }

    /// De-allocate all the memory previously allocated with @alloc
    void freeAll() {
        idx = 0;
    }

    /// Get the max number of bytes that can be allocated by the allocator
    int64_t maxBytes() const {
        return totalBytes;
    }

    /// Get the free space that can still be allocated, same as maxBytes before any allocations
    int64_t freeBytes() const {
        return totalBytes - idx;
    }

    void zeroAll() const {
        memset(data, 0, totalBytes);
    }

    StackAllocator(const StackAllocator &) = delete;
    StackAllocator &operator=(const StackAllocator &) = delete;
private:
    const int64_t totalBytes;
    int64_t idx = 0;
    uint8_t *data = nullptr;
};

void ThrashCache() {
    static std::vector<uint8_t> trash(100_mb);
    std::atomic_thread_fence(std::memory_order_acq_rel);
    volatile uint8_t *data = trash.data();
    for (int c = 0; c < 10; c++) {
        for (int i = 0; i < int(trash.size()); i++) {
            data[i] = data[i] * 2 + 1;
        }
    }
    std::atomic_thread_fence(std::memory_order_acq_rel);
}


/**
 * Interface for the city stats
 * The city stats are used to store information about cities
 * The cities have an id, name, direction, temperature and humidity
 * The direction is one of the following: north, south, east, west, north-east, north-west, south-east, south-west
 * The temperature is a floating point number
 * The humidity is a floating point number
 * The city stats can be loaded from a file, updated with commands and saved to a file
 */
struct CityStatsInterface {
    StackAllocator *allocator = nullptr;
    CityStatsInterface(StackAllocator *allocator) : allocator(allocator) {}
    virtual void LoadFromFile(std::istream &in) = 0;
    virtual void ExecuteCommands(std::istream &commands) = 0;
    virtual void SaveData(std::ostream &temperature, std::ostream &humidity, std::ostream &directions) = 0;
};

const int MAX_CITY_NAME_LEN = 31;

struct CityInfo {
    int64_t id = -1;
    char name[MAX_CITY_NAME_LEN + 1]{};
    char direction[16]{};
    float temperature{};
    float humidity{};
};

float RoundFloat(float value) {
    return float(int(value * 10) / 10.0);
}

/**
 * CityStats implementation
 * Sample implementation of the CityStatsInterface
 */
struct CityStats : CityStatsInterface {
    std::vector<CityInfo> data;
    std::unordered_map<std::string, std::string> leftRotatedDir, rightRotatedDir;

    CityStats(StackAllocator *allocator) : CityStatsInterface(allocator) {
        leftRotatedDir["north"] = "north-west";
        leftRotatedDir["north-west"] = "west";
        leftRotatedDir["west"] = "south-west";
        leftRotatedDir["south-west"] = "south";
        leftRotatedDir["south"] = "south-east";
        leftRotatedDir["south-east"] = "east";
        leftRotatedDir["east"] = "north-east";
        leftRotatedDir["north-east"] = "north";

        rightRotatedDir["north"] = "north-east";
        rightRotatedDir["north-east"] = "east";
        rightRotatedDir["east"] = "south-east";
        rightRotatedDir["south-east"] = "south";
        rightRotatedDir["south"] = "south-west";
        rightRotatedDir["south-west"] = "west";
        rightRotatedDir["west"] = "north-west";
        rightRotatedDir["north-west"] = "north";
    }

    void LoadFromFile(std::istream &in) override {
        CityInfo city;
        while (in >> city.id >> city.name >> city.direction >> city.temperature >> city.humidity) {
            data.push_back(city);
        }
    }

    virtual void ExecuteCommands(std::istream &commands) override {
        char type;
        int64_t from, to;
        float delta;
        char rotate;
        while (commands >> type >> from >> to >> delta >> rotate) {
            if (type == 't') {
                UpdateTemperatureInRange(from, to, delta, rotate);
            } else {
                UpdateHumidityInRange(from, to, delta, rotate);
            }
        }
    }

    void UpdateTemperatureInRange(int64_t startID, int64_t endID, float delta, char rotate) {
        auto start = std::lower_bound(data.begin(), data.end(), startID, [](const CityInfo &city, int64_t id) {
            return city.id < id;
        });
        auto end = std::upper_bound(data.begin(), data.end(), endID, [](int64_t id, const CityInfo &city) {
            return id < city.id;
        });

        for (auto it = start; it != end; ++it) {
            if (rotate == 'r') {
                strcpy(it->direction, rightRotatedDir[it->direction].c_str());
            } else {
                strcpy(it->direction, leftRotatedDir[it->direction].c_str());
            }
            it->temperature += delta;
        }
    }

    void UpdateHumidityInRange(int64_t startID, int64_t endID, float delta, char rotate) {
        auto start = std::lower_bound(data.begin(), data.end(), startID, [](const CityInfo &city, int64_t id) {
            return city.id < id;
        });
        auto end = std::upper_bound(data.begin(), data.end(), endID, [](int64_t id, const CityInfo &city) {
            return id < city.id;
        });

        for (auto it = start; it != end; ++it) {
            if (rotate == 'r') {
                strcpy(it->direction, rightRotatedDir[it->direction].c_str());
            } else {
                strcpy(it->direction, leftRotatedDir[it->direction].c_str());
            }
            it->humidity += delta;
            if (it->humidity < 0) {
                it->humidity = 0;
            } else if (it->humidity > 100) {
                it->humidity = 100;
            }
        }
    }

    virtual void SaveData(std::ostream &temperature, std::ostream &humidity, std::ostream &directions) override {
        float sumTemp = 0;
        float sumHumidity = 0;
        std::unordered_map<std::string, int> dirCount;
        for (const auto &city : data) {
            sumTemp += city.temperature;
            sumHumidity += city.humidity;
            dirCount[city.direction]++;
            temperature.write(reinterpret_cast<const char *>(&city.temperature), sizeof(float));
            humidity.write(reinterpret_cast<const char *>(&city.humidity), sizeof(float));
            directions << city.direction << ' ';
        }
        const float avgHumidity = float(int(sumHumidity / data.size()));
        const float avgTemperature = float(int(sumTemp / data.size()));

        temperature.write(reinterpret_cast<const char *>(&avgTemperature), sizeof(float));
        humidity.write(reinterpret_cast<const char *>(&avgHumidity), sizeof(float));

        std::string mostCommonDir = dirCount.begin()->first;
        for (const auto &dir : dirCount) {
            if (dir.second > dirCount[mostCommonDir]) {
                mostCommonDir = dir.first;
            }
        }
        directions << mostCommonDir << ' ';
    }
};

struct TestDescription {
    int cityCount;
    int64_t itemCount;
    int64_t updateCount;
    std::string name;
    int repeat = 10;
};

struct DataGenerator {
    std::vector<std::string> cities;
    std::vector<std::string> directions = { "north", "south", "east", "west", "north-east", "north-west", "south-east", "south-west" };
    std::mt19937_64 rng{ 0xdeadbeef };

    DataGenerator() {
    }

    void GenerateData(const TestDescription &desc) {
        cities.clear();
        GenerateRandomCities(desc.cityCount);
        int64_t start, end;
        GenerateCityData(desc.name, desc.itemCount, start, end);
        GenerateUpdateCommands(desc.name, desc.updateCount, start, end);
    }

    static void GenerateTestData(const std::vector<TestDescription> &tests) {
        DataGenerator gen;
        for (const auto &desc : tests) {
            std::cout << "Generating data for test " << desc.name << std::endl;
            gen.GenerateData(desc);
        }
    }

    void GenerateUpdateCommands(const std::string &name, int64_t count, int64_t startID, int64_t endID) {
        std::uniform_int_distribution<int64_t> id(startID, endID);
        std::uniform_real_distribution<float> delta(-10, 10);
        std::uniform_int_distribution<int> typePicker(0, 1);
        char rotateDir[2] = { 'r', 'l' };
        char type[2] = { 't', 'h' };

        std::ofstream out(name + "-updates.txt", std::ios::binary);
        for (int64_t c = 0; c < count; c++) {
            const int64_t a = id(rng);
            const int64_t b = id(rng);
            out << type[typePicker(rng)] << ' ' << std::min(a, b) << ' ' << std::max(a, b) << ' ' << RoundFloat(delta(rng)) << ' ' << rotateDir[typePicker(rng)] << '\n';
        }
    }

    void GenerateCityData(const std::string &name, int64_t count, int64_t &startID, int64_t &endID) {
        std::uniform_int_distribution<int> cityPicker(0, int(cities.size() - 1));
        std::uniform_int_distribution<int> idSkip(1, 10);
        std::normal_distribution<float> temperature(20, 10);
        std::uniform_real_distribution<float> humidity(0, 100);
        std::uniform_int_distribution<int> directionPicker(0, int(directions.size() - 1));

        std::ofstream out(name + "-cities.txt", std::ios::binary);
        startID = idSkip(rng);
        for (int64_t c = 0, id = startID; c < count; c++) {
            out << id << ' ' << cities[cityPicker(rng)] << ' ' << directions[directionPicker(rng)] << ' ' << RoundFloat(temperature(rng)) << ' ' << RoundFloat(humidity(rng)) << '\n';
            id += idSkip(rng);
            endID = id;
        }
    }

    void GenerateRandomCities(int count) {
        cities.resize(count);
        std::uniform_int_distribution<int> letter('a', 'z');
        std::uniform_int_distribution<int> cityLen(4, MAX_CITY_NAME_LEN);
        std::generate(cities.begin(), cities.end(), [this, &letter, &cityLen]() {
            std::string city;
            const int len = cityLen(rng);
            for (int i = 0; i < len; ++i) {
                city.push_back(letter(rng));
            }
            return city;
        });
    }
};

struct Tester {
    CityStatsInterface *base;
    CityStatsInterface *update;
    bool timeTester;
    Tester(CityStatsInterface *base, CityStatsInterface *update, bool timeTester) : base(base), update(update), timeTester(timeTester) {}

    struct TestResult {
        bool correct{ false };
        struct Times {
            std::chrono::nanoseconds loadTime{};
            std::chrono::nanoseconds commandsTime{};
            std::chrono::nanoseconds saveTime{};
        };
        Times base;
        Times update;
    };

    TestResult RunTest(const std::string &name) {
        std::vector<float> baseData;
        std::vector<float> updateData;
        std::string baseDir, updateDir;
        TestResult result;

        if (timeTester) {
            ThrashCache();
        }
        result.base = TestImplementation(name, base, baseData, baseDir, "base");
        if (timeTester) {
            ThrashCache();
        }
        result.update = TestImplementation(name, update, updateData, updateDir, "update");

        if (!timeTester) {
            result.correct = FindDiffInData(baseData, updateData) == -1;
            if (baseDir != updateDir) {
                std::cout << "Different directions: [" << baseDir << ']' << std::endl << std::endl << '[' << updateDir << ']' << std::endl;
                result.correct = false;
            }
        }

        return result;
    }

    template <typename Function>
    static std::chrono::nanoseconds TestFunction(Function &&f) {
        auto start = std::chrono::high_resolution_clock::now();
        std::atomic_thread_fence(std::memory_order_acq_rel);
        f();
        std::atomic_thread_fence(std::memory_order_acq_rel);
        return std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - start);
    }

    static TestResult::Times TestImplementation(const std::string &name, CityStatsInterface *test, std::vector<float> &data, std::string &directionsOut, const std::string &type) {
        std::ifstream in(name + "-cities.txt", std::ios::binary);
        std::ifstream commands(name + "-updates.txt", std::ios::binary);
        TestResult::Times res;

        res.loadTime = TestFunction([&]() {
            test->LoadFromFile(in);
        });

        res.commandsTime = TestFunction([&]() {
            test->ExecuteCommands(commands);
        });

        std::stringstream temperature;
        std::stringstream humidity;
        std::stringstream directions;
        res.saveTime = TestFunction([&]() {
            test->SaveData(temperature, humidity, directions);
        });

        directionsOut = directions.str();
        const int64_t tempSize = temperature.tellp() / sizeof(float);
        const int64_t humSize = humidity.tellp() / sizeof(float);
        temperature.seekg(0);
        humidity.seekg(0);

        data.resize(tempSize + humSize);
        temperature.read(reinterpret_cast<char *>(data.data()), tempSize * sizeof(float));
        humidity.read(reinterpret_cast<char *>(data.data() + tempSize), humSize * sizeof(float));

        return res;
    }

    int64_t FindDiffInData(const std::vector<float> &a, const std::vector<float> &b) {
        if (a.size() != b.size()) {
            std::cout << "Difference in size " << a.size() << " " << b.size() << std::endl;
            return -2;
        }
        int64_t diff = -1;
        for (int64_t c = 0; c < int64_t(a.size()); c++) {
            if (a[c] != b[c]) {
                std::cout << "Difference at " << c << " " << a[c] << " " << b[c] << std::endl;
                if (diff == -1) {
                    diff = c;
                }
            }
        }
        return diff;
    }
};

std::ostream &PrintTime(std::ostream &out, const std::chrono::nanoseconds &ns) {
    auto us = std::chrono::duration_cast<std::chrono::microseconds>(ns);
    auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(ns);
    auto s = std::chrono::duration_cast<std::chrono::seconds>(ns);
    if (s > 10s) {
        return out << s;
    } else if (ms > 10ms) {
        return out << ms;
    } else if (us > 10us) {
        return out << us;
    } else {
        return out << ns;
    }
}

std::ostream &operator<<(std::ostream &out, const Tester::TestResult::Times &times) {
    PrintTime(out, times.loadTime) << ' ';
    PrintTime(out, times.commandsTime) << ' ';
    PrintTime(out, times.saveTime);
    return out;
}

const int repeat = 5;
std::vector<TestDescription> tests = {
    {3, 1000, 10, "small", repeat},
    {100, 10'000, 10'000, "medium", repeat},
    {100, 50'000, 20'000, "large", repeat},
    {100, 1'000, 1'000'000, "many-updates", repeat},
    {1000, 10'000'000, 100, "many-records", repeat},
    {100'000, 1'000'000, 100, "many-cities", repeat},
};

//////////////////////// BEGIN IMPLEMENTATION ////////////////////////
// No changes outside of this region are allowed
// Implement the FastCityStats class here
#ifdef _WIN64
#define FORCE_INLINE __forceinline
#define NEVER_INLINE __declspec(noinline)
#else
#define FORCE_INLINE __attribute__((always_inline))
#define NEVER_INLINE __attribute__((noinline))
#endif

#define ENABLE_SIMD
//#define TRACY_ENABLE

#ifdef TRACY_ENABLE
#include "tracy/public/tracy/Tracy.hpp"
#else
#define ZoneScoped
#define ZoneScopedN(...)
#define ZoneValue(...)

#endif
#ifdef ENABLE_SIMD
#include <immintrin.h>
#include "vectorclass.h"
#endif
#include <array>
#include <string_view>
#include <thread>

template <typename T>
struct AlignedArrayPtr {
    void *allocated = nullptr;
    T *aligned = nullptr;
    int capacity = -1;
    int count = 0;


    AlignedArrayPtr() = default;

    AlignedArrayPtr(int count) {
        init(count);
    }

    void init(int newCount) {
        if (newCount > capacity) {
            free(allocated);
            aligned = alignedAlloc<T>(newCount, allocated);
            capacity = newCount;
        }
    }

    void reserve(int newCap) {
        init(newCap);
    }

    void memset(int value) {
        ::memset(aligned, value, sizeof(T) * capacity);
    }

    ~AlignedArrayPtr() {
        free(allocated);
    }

    T *get() {
        return aligned;
    }

    const T *get() const {
        return aligned;
    }

    const T *data() const {
        return get();
    }

    T *data() {
        return get();
    }

    int size() const {
        return count;
    }

    operator T *() {
        return aligned;
    }

    operator const T *() const {
        return aligned;
    }

    int64_t getCount() const {
        return capacity;
    }

    const T *begin() const {
        return aligned;
    }

    const T *end() const {
        return aligned + count;
    }

    const T &operator[](int index) const {
        return aligned[index];
    }

    T &operator[](int index) {
        return aligned[index];
    }

    const T &operator[](int64_t index) const {
        return aligned[index];
    }

    T &operator[](int64_t index) {
        return aligned[index];
    }

    void push_back(const T &item) {
        aligned[count++] = item;
    }

    void resize(int newCount) {
        count = newCount;
    }

    void clear() {
        count = 0;
    }

    AlignedArrayPtr(const AlignedArrayPtr &) = delete;
    AlignedArrayPtr &operator=(const AlignedArrayPtr &) = delete;
};

[[noreturn]] inline void unreachable() {
    // Uses compiler specific extensions if possible.
    // Even if no extension is used, undefined behavior is still raised by
    // an empty function body and the noreturn attribute.
#if defined(_MSC_VER) && !defined(__clang__) // MSVC
    __assume(false);
#else // GCC, Clang
    __builtin_unreachable();
#endif
}

struct FastCityStats : CityStatsInterface {
    enum Directions : uint8_t {
        NORTH,
        NORTH_EAST,
        EAST,
        SOUTH_EAST,
        SOUTH,
        SOUTH_WEST,
        WEST,
        NORTH_WEST,
        _Count,
    };
    struct Command {
        int start;
        bool isRightRotate;
        int end;
        float delta;
    };

    inline static constexpr int RESERVE_SIZE = 10'000'000;

    struct CommandList {
        AlignedArrayPtr<int> start;
        AlignedArrayPtr<int> end;
        AlignedArrayPtr<float> delta;
        AlignedArrayPtr<bool> isRightRotate;

        void reserve() {
            start.reserve(RESERVE_SIZE);
            end.reserve(RESERVE_SIZE);
            delta.reserve(RESERVE_SIZE);
            isRightRotate.reserve(RESERVE_SIZE);
        }

        void clear() {
            start.resize(0);
            end.resize(0);
            delta.resize(0);
            isRightRotate.resize(0);
        }

        int size() const {
            return start.size();
        }

        void push_back(const Command &cmd) {
            ZoneScoped;
            start.push_back(cmd.start);
            end.push_back(cmd.end);
            delta.push_back(cmd.delta);
            isRightRotate.push_back(cmd.isRightRotate);
        }
    };

    inline static AlignedArrayPtr<int> dataID;
    inline static AlignedArrayPtr<float> dataTemp;
    inline static AlignedArrayPtr<float> dataHumidity;
    inline static AlignedArrayPtr<uint8_t> dataDirs;
    inline static AlignedArrayPtr<char> inputData;
    inline static CommandList hCommands, tCommands;

    char dirNamesBuff[68] = { "north north-east east south-east south south-west west north-west " };
    std::string_view dirNames[8] = {
        {dirNamesBuff + 0, 6},
        {dirNamesBuff + 6, 11},
        {dirNamesBuff + 17, 5},
        {dirNamesBuff + 22, 11},
        {dirNamesBuff + 33, 6},
        {dirNamesBuff + 39, 11},
        {dirNamesBuff + 50, 5},
        {dirNamesBuff + 55, 11},
    };
    inline static constexpr int READ_BUFFER_SIZE = 450_mb;
    inline static constexpr int ID_MAX_DIGIT = 21;

    FastCityStats(StackAllocator *allocator) : CityStatsInterface(allocator) {
        ZoneScoped;
        dataID.reserve(RESERVE_SIZE);
        dataTemp.reserve(RESERVE_SIZE);
        dataHumidity.reserve(RESERVE_SIZE);
        dataDirs.reserve(RESERVE_SIZE);
        inputData.reserve(READ_BUFFER_SIZE + 10);
        hCommands.reserve();
        tCommands.reserve();

        hCommands.clear();
        tCommands.clear();

        dataID.resize(0);
        dataTemp.resize(0);
        dataHumidity.resize(0);
        dataDirs.resize(0);
        inputData.resize(0);
    }

    FORCE_INLINE static int ParseInt(const char *str, int maxlen, int &index) {
        ZoneScoped;
        int value = 0;
        const int start = index;
        for (int c = 0; c < maxlen && str[index] != ' ' && str[index] != '\n' && str[index] != '.'; index++, c++) {
            value = value * 10 + str[index] - '0';
        }
        ZoneValue(value);
        return value;
    }

    FORCE_INLINE static int ParseIntSimd(const char *str, int maxlen, int &index) {
        ZoneScoped;
#ifdef ENABLE_SIMD
        constexpr int multipliers[16] = {
            10'000'000,
            1'000'000,
            100'000,
            10'000,
            1'000,
            100,
            10,
            1,
        };

        const Vec16c maskSpace(' '), maskNewLine('\n'), maskDot('.');
        Vec16c digits;
        digits.load(str + index);

        Vec16cb stopMask = (digits == maskSpace) | (digits == maskNewLine) | (digits == maskDot);
        const int stopIndex = BSF(uint32_t(to_bits(stopMask)));
        digits = digits - Vec16c('0');
        Vec8i powers, digitNumbers = _mm256_cvtepi8_epi32(digits);
        powers.load(multipliers + (8 - stopIndex));

        auto parts = powers * digitNumbers;
        index += stopIndex;
        const int value = horizontal_add_x(parts);
        ZoneValue(value);
        return value;

#else
        return ParseInt(str, maxlen, index);
#endif
    }


    FORCE_INLINE static Directions GetDir(const char *dirName) {
        ZoneScoped;
 
        int k1 = dirName[0];
        int k2 = dirName[5] == ' ';
        int k3 = dirName[6] == 'w';
        int k4 = dirName[6] == 'e';

        int res
            = (k1 == 'e') * EAST
            + (k1 == 'w') * WEST
            + (k1 == 'n') * k2 * NORTH
            + (k1 == 'n') * k3 * NORTH_WEST
            + (k1 == 'n') * k4 * NORTH_EAST
            + (k1 == 's') * k2 * SOUTH
            + (k1 == 's') * k3 * SOUTH_WEST
            + (k1 == 's') * k4 * SOUTH_EAST;

        return Directions(res);
    }

    template <char stop>
    FORCE_INLINE float ParseFloat(const char *str, int maxLen, int &index) {
        ZoneScoped;
        const int negative = str[index] == '-';
        index += negative;
        int64_t nominator = ParseInt(str, maxLen, index);
        const float multiplier[] = { 1.f, -1.f };
        const bool hasDenom = str[index] == '.';
        const uint8_t denomDigit = str[index + 1] - '0';
        const float denomSelector[2] = { 0.f, 1.f };
        const float denomMask = denomSelector[hasDenom];
        const float denom = (denomDigit / 10.f) * denomMask;
        const float value = (nominator + denom) * multiplier[negative];
        index += int(denomMask) * 2;

        ZoneValue(value);
        return value;
    }

    FORCE_INLINE static int BSF(uint64_t value) {
#ifdef _WIN64
        unsigned long int res;
        _BitScanForward64(&res, value);
        return res;
#else
        return __builtin_ctzll(value);
#endif
    }

    FORCE_INLINE static int BSF(uint32_t value) {
#ifdef _WIN64
        unsigned long int res;
        _BitScanForward(&res, value);
        return res;
#else
        return __builtin_ctz(value);
#endif
    }

    FORCE_INLINE void ParseLine(int start, int length, int &longIDs) {
        ZoneScoped;
        std::string_view line(inputData.data() + start, length);

        const char *lineData = inputData.data() + start;

        int index = 0;
        int id;
        if (longIDs) {
            id = ParseIntSimd(lineData, ID_MAX_DIGIT, index);
        } else {
            id = ParseInt(lineData, ID_MAX_DIGIT, index);
            longIDs = index > 4;
        }
        index++;
#ifdef ENABLE_SIMD
        {
            ZoneScopedN("SkipCityName");
            Vec32c inputData, spaceMask(' ');
            inputData.load(lineData + index);
            const uint32_t mask = to_bits(inputData == spaceMask);
            const int spaceIndex = BSF(mask);
            index = index + spaceIndex + 1;
        }

#else
        for (int c = 0; c < 32 && lineData[index] != ' '; c++) {
            index++;
        }
        index++;
#endif

        const auto dir = GetDir(lineData + index);
        {
            ZoneScopedN("SkipDirName");
            Vec32c inputData, spaceMask(' ');
            inputData.load(lineData + index);
            const uint32_t mask = to_bits(inputData == spaceMask);
            const int spaceIndex = BSF(mask);
            index = index + spaceIndex + 1;
        }
        //index += dirNames[dir].length();

        const float temp = ParseFloat<' '>(lineData, 2, index);
        ++index;
        const float hum = ParseFloat<'\n'>(lineData, 2, index);

        dataID.push_back(id);
        dataTemp.push_back(temp);
        dataHumidity.push_back(hum);
        dataDirs.push_back(dir);
    }

    void LoadFromFile(std::istream &in) override {
        ZoneScoped;
        int longIds = false;
        ParseInputLines<3, 8>(in, [&](int start, int length) {
            ParseLine(start, length, longIds);
        });
    }

    FORCE_INLINE void ParseLineCommand(int64_t start, int64_t length) {
        ZoneScoped;
        const char *lineData = inputData.data() + start;
        std::string_view line(lineData, length);
        int index = 0;
        const char type = lineData[index];
        index += 2;
        Command cmd;
        cmd.start = ParseInt(lineData, ID_MAX_DIGIT, index);
        index++;
        cmd.end = ParseInt(lineData, ID_MAX_DIGIT, index);
        index++;
        cmd.delta = ParseFloat<' '>(lineData, 2, index);
        index++;
        const char rot = lineData[index++];
        cmd.isRightRotate = rot == 'r';
        if (type == 't') {
            tCommands.push_back(cmd);
        } else {
            hCommands.push_back(cmd);
        }
    }

    inline static constexpr int INDEX_SIZE = 1024;
    using IndexContainer = std::array<int, INDEX_SIZE + 2>;

    FORCE_INLINE void buildIDIndex(IndexContainer &index) {
        ZoneScoped;
        const float elementPerInterval = float(dataID.size()) / INDEX_SIZE;
        for (int c = 0; c < INDEX_SIZE; c++) {
            index[c + 1] = dataID[std::min<int>((c + 1) * elementPerInterval, dataID.size() - 1)];
        }
        index[0] = dataID[0];
        index[index.size() - 1] = dataID[dataID.size() - 1];
    }

    bool isInRange(int a, int b, int val) {
        return val >= a && val <= b;
    }

    FORCE_INLINE void remapIndices(CommandList &commands, IndexContainer &index) {
        ZoneScoped;
        const float elementPerInterval = float(dataID.size()) / INDEX_SIZE;

        const int len = commands.size();
        const int idCount = dataID.size();
        for (int c = 0; c < len; c++) {
            const int startID = commands.start[c];
            const int endID = commands.end[c];

            const int startIndexIndex = std::lower_bound(index.begin(), index.end(), startID) - index.begin();
            const int endIndexIndex = std::upper_bound(index.begin(), index.end(), endID) - index.begin();

            //assert(isInRange(index[std::clamp<int>(startIndexIndex - 1, 0, index.size() - 1)], index[std::clamp<int>(startIndexIndex + 1, 0, index.size() - 1)], startID)
            //    || startID > index.back() || startID < index.front());
            //assert(isInRange(index[std::clamp<int>(endIndexIndex - 1, 0, index.size() - 1)], index[std::clamp<int>(endIndexIndex + 1, 0, index.size() - 1)], endID)
            //    || endID > index.back() || endID < index.front());

            auto start = std::lower_bound(
                dataID.begin() + std::max<int>(0, (startIndexIndex - 2) * elementPerInterval),
                dataID.begin() + std::min<int>(idCount, (startIndexIndex + 0) * elementPerInterval),
                startID);
            //auto _start = std::lower_bound(dataID.begin(), dataID.end(), startID);
            //assert(start == _start);

            auto end = std::upper_bound(
                dataID.begin() + std::max<int>(0, (endIndexIndex - 2) * elementPerInterval),
                dataID.begin() + std::min<int>(idCount, (endIndexIndex + 0) * elementPerInterval),
                endID);
            //auto _end = std::upper_bound(dataID.begin(), dataID.end(), endID);
            //assert(end == _end);

            auto startIndex = start - dataID.begin();
            auto endIndex = end - dataID.begin();
            commands.start[c] = startIndex;
            commands.end[c] = endIndex;
        }
    }

    template <int loopCount = 8, int arrSize>
    static FORCE_INLINE int FindLineEndings(char *data, std::array<int, arrSize> &indices) {
        ZoneScoped;
        int count = 0;
#ifdef ENABLE_SIMD
        Vec64c mask('\n');
        Vec64c dataVec;
        for (int c = 0; c < loopCount; c++) {
            dataVec.load(data);
            data += 64;
            const Vec64c result = dataVec == mask;
            uint64_t maskResult = to_bits(result);
            const int bits = vml_popcnt(maskResult);
            for (int r = 0; r < bits; r++) {
                const int pos = BSF(maskResult);
                indices[count++] = (c * 64) + (pos + 1);
                maskResult ^= (1llu << pos);
            }
        }
        ZoneValue(count);
#endif
        return count;
    }

    template <int MaxLinePer64Byte = 4, int LoopCount = 16, typename ParseFN>
    static FORCE_INLINE void ParseInputLines(std::istream &in, ParseFN &&fn) {
        ZoneScoped;
        int64_t readSize;
        {
            ZoneScopedN("Read from input stream");
            in.read(inputData.data(), READ_BUFFER_SIZE);
            readSize = in.gcount();
            inputData[readSize] = '\n';
        }

        for (int q = 0; q < 100; q++) {
            int prev = 0;
            int batchCount = 0;
            const int loopCount = LoopCount;
#ifdef ENABLE_SIMD
            const int arrSize = loopCount * MaxLinePer64Byte;
            batchCount = readSize / (loopCount * 64);
            std::array<int, arrSize> indices{};
            for (int b = 0; b < batchCount; b++) {
                const int offset = b * (loopCount * 64);
                const int count = FindLineEndings<loopCount, arrSize>(inputData.data() + offset, indices);
                for (int c = 0; c < count; c++) {
                    const int start = prev;
                    const int end = indices[c];
                    prev = end;
                    fn(start + offset, end - start);
                }
                prev = -((loopCount * 64) - prev);
            }
#endif
            int64_t index = batchCount * (loopCount * 64) + prev;

            while (true) {
                const int64_t start = index;
                while (index < readSize && inputData[index] != '\n') {
                    ++index;
                }

                if (index - start == 0) {
                    break;
                }
                fn(start, index - start);
                ++index;
            }
            return;
            dataID.clear();
            dataTemp.clear();
            dataHumidity.clear();
            dataDirs.clear();
        }
    }

    virtual void ExecuteCommands(std::istream &commands) override {
        ZoneScoped;
        ParseInputLines(commands, [this](int start, int length) {
            ParseLineCommand(start, length);
        });

        {
            ZoneScopedN("RemapIndices");
            IndexContainer index;
            buildIDIndex(index);
            remapIndices(hCommands, index);
            remapIndices(tCommands, index);
        }


        char dir[2]{ 'l', 'r' };
        const int humCmdCount = int(hCommands.size());
        const int tempCmdCount = int(tCommands.size());
        for (int c = 0; c < humCmdCount; c++) {
            UpdateHumidityInRange(hCommands.start[c], hCommands.end[c], hCommands.delta[c], dir[hCommands.isRightRotate[c]]);
        }
        for (int c = 0; c < tempCmdCount; c++) {
            UpdateTemperatureInRange(tCommands.start[c], tCommands.end[c], tCommands.delta[c], dir[tCommands.isRightRotate[c]]);
        }
    }

    FORCE_INLINE void ChangeDirs(int64_t startIndex, int64_t endIndex, int8_t change) {
        ZoneScoped;
        ZoneValue(endIndex - startIndex);
        const uint8_t maskValue = _Count - 1;
        int64_t c = startIndex;
#ifdef ENABLE_SIMD
        Vec32c dirData;
        const Vec32c changeValue(change);

        const Vec32c mask(maskValue);
        const int untilAligned = std::min<int>(endIndex - startIndex, (32 - (intptr_t(dataDirs.data() + c) % 32)) / sizeof(uint8_t));

        for (int r = 0; r < untilAligned; r++) {
            dataDirs[r + c] = uint8_t(dataDirs[r + c] + change) & maskValue;
        }
        c += untilAligned;

        while (endIndex - c >= 32) {
            dirData.load(dataDirs.data() + c);
            dirData = (dirData + changeValue) & mask;
            dirData.store(dataDirs.data() + c);
            c += 32;
        }
#endif

        for (; c < endIndex; c++) {
            dataDirs[c] = uint8_t(dataDirs[c] + change) & maskValue;
        }
    }

    FORCE_INLINE void UpdateTemperatureInRange(int64_t startID, int64_t endID, float delta, char rotate) {
        ZoneScoped;
        ZoneValue(endID - startID);
        auto startIndex = startID;
        auto endIndex = endID;
        int64_t c = startIndex;
#ifdef ENABLE_SIMD

        Vec8f deltaV = delta, data;
        const int untilAligned = std::min<int>(endIndex - startIndex, (32 - (intptr_t(dataTemp.data() + c) % 32)) / sizeof(float));
        for (int r = 0; r < untilAligned; r++) {
            dataTemp[r + c] += delta;
        }
        c += untilAligned;

        while (endIndex - c >= 8) {
            data.load_a(dataTemp.data() + c);
            data += deltaV;
            data.store_a(dataTemp.data() + c); // store_nt
            c += 8;
        }
#endif

        for (; c < endIndex; c++) {
            dataTemp[c] += delta;
        }

        ChangeDirs(startIndex, endIndex, rotate == 'r' ? 1 : -1);
    }

    FORCE_INLINE void UpdateHumidityInRange(int64_t startID, int64_t endID, float delta, char rotate) {
        ZoneScoped;
        ZoneValue(endID - startID);
        auto startIndex = startID;
        auto endIndex = endID;

        int64_t c = startIndex;
#ifdef ENABLE_SIMD
        Vec8f deltaV = delta, data;
        Vec8f maxHum = 100, minHum = 0;
        const int untilAligned = std::min<int>(endIndex - startIndex, (32 - (intptr_t(dataHumidity.data() + c) % 32)) / sizeof(float));
        for (int r = 0; r < untilAligned; r++) {
            dataHumidity[r + c] = std::clamp(dataHumidity[r + c] + delta, 0.f, 100.f);;
        }
        c += untilAligned;

        while (endIndex - c >= 8) {
            data.load_a(dataHumidity.data() + c);

            data += deltaV;
            data = max(min(data, maxHum), minHum);

            data.store_a(dataHumidity.data() + c); // store_nt
            c += 8;
        }
#endif
        for (; c < endIndex; c++) {
            dataHumidity[c] = std::clamp(dataHumidity[c] + delta, 0.f, 100.f);
        }

        ChangeDirs(startIndex, endIndex, rotate == 'r' ? 1 : -1);
    }

    virtual void SaveData(std::ostream &temperature, std::ostream &humidity, std::ostream &directions) override {
        ZoneScoped;
        temperature.write(reinterpret_cast<const char *>(dataTemp.data()), dataTemp.size() * sizeof(float));
        humidity.write(reinterpret_cast<const char *>(dataHumidity.data()), dataHumidity.size() * sizeof(float));
        {
            ZoneScopedN("ComputeAverages");
            float totalTemp = 0;
            float totalHumidity = 0;
            int64_t c = 0;
#ifdef ENABLE_SIMD
            Vec8f tempAcc = 0, humAcc = 0;
            Vec8f curTemp, curHum;

            for (c = 0; c < int(dataTemp.size()); c += 8) {
                curTemp.load_a(dataTemp.data() + c);
                curHum.load_a(dataHumidity.data() + c);

                tempAcc += curTemp;
                humAcc += curHum;
            }

            for (int r = 0; r < 8; r++) {
                totalTemp += tempAcc[r];
                totalHumidity += humAcc[r];
            }
#endif
            for (; c < dataTemp.size(); c++) {
                totalTemp += dataTemp[c];
                totalHumidity += dataHumidity[c];
            }

            const float avgTemp = int(totalTemp / dataTemp.size());
            const float avgHumidity = int(totalHumidity / dataHumidity.size());
            temperature.write(reinterpret_cast<const char *>(&avgTemp), sizeof(float));
            humidity.write(reinterpret_cast<const char *>(&avgHumidity), sizeof(float));
        }

        int64_t dirCount[_Count] = {};
        {
            ZoneScopedN("Write dirs");
            for (int c = 0; c < int(dataDirs.size()); c++) {
                dirCount[dataDirs[c]]++;
                std::string_view dir = dirNames[dataDirs[c]];
                directions.write(dir.data(), dir.size());
            }
        }


        int mostCommonDir = 0;
        for (int c = 1; c < _Count; c++) {
            if (dirCount[c] > dirCount[mostCommonDir]) {
                mostCommonDir = c;
            }
        }
        directions.write(dirNames[mostCommonDir].data(), dirNames[mostCommonDir].size());
    }
};
//////////////////////// END IMPLEMENTATION ////////////////////////

int main(int argc, char *argv[]) {
    StackAllocator empty(nullptr, 0);
    StackAllocator allocator(new uint8_t[1_gb], 1_gb);

    if (0) {
        allocator.freeAll();
        FastCityStats update(&allocator);
        std::vector<float> data;
        std::string output;
        Tester::TestImplementation(tests[4].name, &update, data, output, "update");
        return 0;
    }

    if (argc == 2) {
        if (strcmp(argv[1], "generate") == 0) {
            DataGenerator::GenerateTestData(tests);
        } else if (strcmp(argv[1], "check") == 0) {
            for (const auto &desc : tests) {
                CityStats base(&empty);
                FastCityStats update(&allocator);
                Tester test(&base, &update, false);
                std::cout << "Results for: \"" << desc.name << "\":" << std::endl;
                const auto res = test.RunTest(desc.name);
                if (res.correct) {
                    std::cout << "- Correct" << std::endl;
                } else {
                    std::cout << "- Incorrect" << std::endl;
                }
            }
        } else if (std::isdigit(argv[1][0])) {
            const int idx = atoi(argv[1]);
            if (idx < 0 || idx >= int(tests.size())) {
                std::cout << "Invalid test index\n";
                return 1;
            }
            std::vector<Tester::TestResult::Times> res;
            res.reserve(repeat);
            std::vector<float> data;
            std::string output;
            std::cout << "Running test " << tests[idx].name << " " << repeat << " times\n";
            for (int c = 0; c < repeat; c++) {
                allocator.freeAll();
                FastCityStats update(&allocator);
                const auto times = Tester::TestImplementation(tests[idx].name, &update, data, output, "update");
                res.push_back(times);
            }
            for (const auto &time : res) {
                std::cout << time << std::endl;
            }
        } else if (!strcmp(argv[1], "help")) {
            std::cout << "Usage: CityInfo [command]\n";
            std::cout << "Commands:\n";
            std::cout << "  generate - generate test data\n";
            std::cout << "  check - check the correctness of the implementation\n";
            std::cout << "  [index] - run a single test\n";
            std::cout << "  <no-args> - run all timing tests\n";
        } else {
            std::cout << "Unknown command\n";
        }
        std::cout << "Waiting for 5 seconds before exiting" << std::endl;
        std::this_thread::sleep_for(1s);
        return 0;
    }

    std::shuffle(tests.begin(), tests.end(), std::mt19937_64(std::random_device()()));

    for (const auto &desc : tests) {
        std::cout << "Results for " << desc.name << std::endl;
        std::cout << "base(load, commands, save) | update(load, commands, save)" << std::endl;
        const int repeatCount = desc.repeat;
        std::vector<Tester::TestResult::Times> baseResults, updateResults;
        for (int c = 0; c < repeatCount; c++) {
            allocator.zeroAll();
            allocator.freeAll();
            CityStats base(&empty);
            FastCityStats update(&allocator);
            Tester test(&base, &update, true);
            const auto res = test.RunTest(desc.name);
            baseResults.push_back(res.base);
            updateResults.push_back(res.update);
            std::cout << baseResults[c] << " | ";
            std::cout << updateResults[c] << std::endl;
        }
        int bestBase = 0, bestUpdate = 0;
        std::chrono::nanoseconds totalBase{}, totalUpdate{};
        for (int c = 0; c < int(baseResults.size()); c++) {
            const auto baseTime = baseResults[c].loadTime + baseResults[c].commandsTime + baseResults[c].saveTime;
            const auto updateTime = updateResults[c].loadTime + updateResults[c].commandsTime + updateResults[c].saveTime;
            totalBase += baseTime;
            totalUpdate += updateTime;
            if (baseTime < baseResults[bestBase].loadTime + baseResults[bestBase].commandsTime + baseResults[bestBase].saveTime) {
                bestBase = c;
            }
            if (updateTime < updateResults[bestUpdate].loadTime + updateResults[bestUpdate].commandsTime + updateResults[bestUpdate].saveTime) {
                bestUpdate = c;
            }
        }
        std::cout << "Average base: ";
        PrintTime(std::cout, totalBase / repeatCount) << " Average update: ";
        PrintTime(std::cout, totalUpdate / repeatCount) << std::endl;
        std::cout << "Best base: " << baseResults[bestBase] << " Best update: " << updateResults[bestUpdate] << std::endl;
        std::cout << std::endl;
        std::cout << "Best load speedup: " << float(baseResults[bestBase].loadTime.count()) / updateResults[bestUpdate].loadTime.count() << "x\n";
        std::cout << "Best commands speedup: " << float(baseResults[bestBase].commandsTime.count()) / updateResults[bestUpdate].commandsTime.count() << "x\n";
        std::cout << "Best save speedup: " << float(baseResults[bestBase].saveTime.count()) / updateResults[bestUpdate].saveTime.count() << "x\n";
        std::cout << std::endl;
        const float speedup = float(totalBase.count()) / totalUpdate.count();
        std::cout << "Total Speedup: " << speedup << "x\n";
        std::cout << "----------------------------------------\n";
    }

    return 0;
}