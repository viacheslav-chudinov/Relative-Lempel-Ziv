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
using int128_t = __int128;


constexpr uint_t NONE = std::numeric_limits<uint_t>::max();

inline constexpr uint128_t SUPER_PRIME = (uint128_t(1)<<127)-1;  //It is the 12th known Mersenne prime
inline constexpr u32 SIGMA = 4;
inline char ALPHABET[SIGMA] = {'a','b','c','d'};


inline std::array<int, 0xff> char_to_index;

inline constexpr uint128_t q_prime = SUPER_PRIME;

// a mod q :
// uint128_t reduced_mod = ((a>>127)&1 + (a&q))
// reduced_mod = reduced_mod - (q&(static_cast<int128_t>(reduced_mod>=q)<<127)>>127) 

// a*b mod q : a = 2^64*d_a + r_a; b = 2^64*d_b + r_b => a*b = 2*d_a*d_b + 2^64 * (d_a*r_b + d_b*r_a) + r_a*r_b
// d_a*r_b + d_b*r_a = 2^64*w + r; => 2^64 * (d_a*r_b + d_b*r_a) = 2*w + 2^64 * r mod q 
//

inline uint128_t reduce(uint128_t a){
    uint128_t t = ((a >> 127) & 1) + (a & q_prime);
    return t - (q_prime & -(uint128_t)(t >= q_prime));
}

inline uint128_t addmod(uint128_t a, uint128_t b){
    return reduce(reduce(a)+reduce(b));
}

inline uint128_t mulmod(uint128_t a, uint128_t b){
    uint128_t d_a = (a >>64) & (~(uint64_t(0)));
    uint128_t d_b = (b >>64) & (~(uint64_t(0)));
    uint128_t r_a = a & (~(uint64_t(0)));
    uint128_t r_b = b & (~(uint64_t(0)));
    uint128_t res = reduce(reduce((reduce(d_a*d_b))<<1)+reduce(r_a*r_b));
    uint128_t intermediate_res = reduce(reduce(d_a*r_b) + reduce(d_b*r_a));
    uint128_t hi = (intermediate_res >>64) & (~(uint64_t(0)));
    uint128_t lo = intermediate_res & (~(uint64_t(0)));
    intermediate_res = reduce((lo<< 64));
    intermediate_res = reduce(intermediate_res + (hi<<1));
    return reduce(intermediate_res + res);
}

inline uint128_t powmod(uint128_t a, uint32_t pow){ // Compute a^pow mod q_prime
    uint128_t tmp = 1;
    uint128_t base = reduce(a);
    while (pow){
        if (pow & 1) tmp = mulmod(tmp, base);
        base = mulmod(base, base);
        pow >>= 1;
    }
    return tmp;
}


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