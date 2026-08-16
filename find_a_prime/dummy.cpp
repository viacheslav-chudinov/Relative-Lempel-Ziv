#include <iostream> 
#include <cmath>
#include "../common_utils.hpp"
#include <iomanip>
#include <bitset>
#include "primes.hpp"

int main(){
    uint32_t m = 9;
    std::cout << std::setprecision(15) << std::log2(m) << std::endl;
    uint128_t limit = from_string("3317044064679887385961981");
    std::cout << limit << std::endl;
    uint32_t n = 1;
    int leading_zeros = __builtin_clz(n);
    int msb_index = 32 - leading_zeros;
    std::cout << "msb_index = " << msb_index << std::endl;
    m = 1<<31;
    uint128_t tmp = static_cast<uint128_t>(m);
    uint128_t k = ((tmp*tmp)*(tmp*tmp)) + 1;
    std::cout << k << std::endl;
    uint64_t j;
    uint32_t r = 1<<30;
    j = (r<<4);
    std::cout << "j = " << std::bitset<64>(j) << std::endl;
    uint128_t p = (uint128_t(1) << 100) - 1;
    for (int i = 0; i < 30; ++i){
        std::cout << "Random in range large p : " << random_in_range_prime(p) << std::endl;
    }
    uint32_t N = 1364723;
    uint64_t q = find_good_prime_for_rb_tables(N);
    std::cout << "A good priime for N = " << N << " is q = "<< q << std::endl;
    std::cout << "q%N =  " << q%N << std::endl;
    uint32_t M = 1776102;
    auto [P, base] =  find_good_prime_for_KR(M, 2);
    std::cout << "A good priime for M is " << P << " and the base is = "<< base << std::endl;
}
