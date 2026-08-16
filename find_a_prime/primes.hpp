#pragma once
#include "../common_utils.hpp"
#include <unordered_set>
#include <cmath>

// if n < 2^64, it is enough to test a = 2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, and 37.
bool is_prime(uint128_t n);

std::pair<uint128_t,uint128_t> find_good_prime_for_KR(uint32_t m, uint32_t lambda);

uint64_t find_good_prime_for_rb_tables(uint32_t N);

uint32_t find_prime_Bertrand(uint32_t i); //Finds a prime in range [i,2i)

