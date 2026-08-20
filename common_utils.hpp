#pragma once
#include <array>
#include <random>
#include <cstdint>
#include <limits>
#include <vector>
#include <iostream>
#include <boost/multiprecision/cpp_int.hpp>

typedef uint32_t u32;
typedef uint64_t u64;
using uint_t = uint64_t;
using uint128_t = unsigned __int128;

constexpr uint_t NONE = std::numeric_limits<uint_t>::max();

inline constexpr uint128_t SUPER_PRIME = (uint128_t(1)<<127)-1;  //It is the 12th known Mersenne prime
inline constexpr u32 SIGMA = 4;
inline char ALPHABET[SIGMA] = {'a','b','c','d'};


inline std::array<int, 0xff> char_to_index;


// NECESSARY TO RUN IN MAIN() BEFORE EVERYTHING:
inline void init_char_to_index(){
    char_to_index.fill(-1);
    for (uint32_t i = 0; i < SIGMA; ++i){
        char_to_index[(unsigned char)ALPHABET[i]] = i;
    }
}

inline uint128_t random_uint128_raw(){
    static std::random_device rd;
    static std::mt19937_64 gen(rd());
    uint128_t high = gen();
    uint128_t low = gen();
    return (high << 64) | low;
}

inline uint128_t random_in_range_prime_with_0(uint128_t prime){
    uint128_t limit = (std::numeric_limits<uint128_t>::max() / prime) * prime; // Introduced to evenly map a large range onto [0..p-1]
    uint128_t x;
    do {
        x = random_uint128_raw();
    } while (x >= limit);
    return x % prime;
}

inline uint128_t random_in_range_prime(uint128_t prime){
    uint128_t range = prime - 1;
    uint128_t limit = (std::numeric_limits<uint128_t>::max() / range) * range;
    uint128_t x;
    do {
        x = random_uint128_raw();
    } while (x >= limit);
    return 1 + (x % range);
}

// inline uint128_t random_in_range_prime(uint128_t prime){
//     static std::random_device rd;
//     static std::mt19937_64 gen(rd());
//     std::uniform_int_distribution<uint128_t> dist(1, prime - 1);
//     return dist(gen);
// }

// inline uint128_t random_in_range_prime_with_0(uint128_t prime){
//     static std::random_device rd;
//     static std::mt19937_64 gen(rd());
//     std::uniform_int_distribution<uint128_t> dist(0, prime - 1);
//     return dist(gen);
// }

inline uint32_t random_in_range(uint32_t N){
    static std::random_device rd;
    static std::mt19937_64 gen(rd());
    std::uniform_int_distribution<uint32_t> dist(0, N - 1);
    return dist(gen);
}

inline uint64_t compute_pow_modulo_n(uint32_t m, uint32_t pow, uint64_t n){
// Compute m^pow modulo n
    uint128_t tmp = 1;
    uint128_t base = ((uint128_t)m) % n;
    while (pow){
        if (pow & 1) {
            tmp = (tmp * base) % n;
        }
        base = (base * base) % n;
        pow >>= 1;
    }
    return (uint64_t)tmp;
}

inline uint128_t compute_pow_modulo_n_long(uint128_t m, uint32_t pow, uint128_t n){
// Compute m^pow modulo n
    using boost::multiprecision::uint256_t;
    uint256_t tmp = 1;
    uint256_t base = m % n;
    uint256_t N = n;
    while (pow){
        if (pow & 1) tmp = (tmp * base) % N;
        base = (base * base) % N;
        pow >>= 1;
    }
    return static_cast<uint128_t>(tmp);
}

inline uint128_t from_string(const std::string& s) {
    uint128_t result = 0;
    for (char c : s) result = result * 10 + (c - '0');
    return result;
}

inline std::ostream& operator<<(std::ostream& os, uint128_t value){
    if (value == 0){
        return os << "0";
    }

    char buf[40];
    char* p = buf + sizeof(buf) - 1;
    while (value != 0){
        *p = static_cast<int>(value%10) + '0';
        p--;
        value = (value-value%10)/10;
    }
    return os.write(p + 1, buf + sizeof(buf) - 1 - p);
}