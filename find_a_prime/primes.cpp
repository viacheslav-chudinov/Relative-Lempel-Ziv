#include "primes.hpp"

// The number of unsuitable bases for Karp-Rabin hash is bounded by 1/24 * m^4 for every m >= 3
// In order to obtain the probability of success of at least 0.75, the length of R cannot be greater than 1_776_102.
bool is_prime(uint128_t n){
    if (!(n&1)){
        return false;
    }
    if (n == 1){
        return false;
    }
    if (n == 2 || n == 3){
        return true;
    }
    std::unordered_set<uint32_t> bases = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41};
    uint128_t tmp = n - 1;
    uint128_t d = 0;
    int s;
    for (s = 0; s < 128; s++){
        if (tmp&1){
            d = tmp;
            break;
        }
        tmp >>= 1;
    }
    if ((n - 1) != (d<<s)){
        throw std::runtime_error("d and s are wrong");
    }
    bool ans = true;
    for (auto a : bases){
        if (a == n){
            return true;
        }
        auto a_pow_d = compute_pow_modulo_n_long(a, d, n);
        if (a_pow_d == 1){
            continue;
        }
        uint128_t tmp = a_pow_d;
        bool found = false;
        for (int r = 0; r < s; r++){
            if (tmp == (n - 1)){
                found = true;
                break;
            }
            tmp = (tmp * tmp)%n;
        }
        if (!found){
            ans = false;
            break;
        }
        
    }
    return ans;
}


std::pair<uint128_t,uint128_t> find_good_prime_for_KR(uint32_t m, uint32_t lambda){  // m is the length of the reference R.
    if(lambda >= 64){
        throw std::runtime_error("Lamda is too large. Execution aborted. It should be that lambda < 64");
    }
    if (m <= 2){
        throw std::runtime_error("R is too short, minimal legit length is 3");
    }
    uint128_t limit = from_string("79609057552317297263087520");
    uint32_t leading_zeros = __builtin_clz(m);
    uint32_t msb_index = 32 - leading_zeros;
    if (4*msb_index + lambda + 1 > 120){
        throw std::runtime_error("Lamda is too large. Execution aborted");
    }
    uint128_t tmp = static_cast<uint128_t>(m);
    uint128_t tmp2 = tmp*tmp;
    uint128_t tmp3 = (((tmp2*tmp2)<<lambda)>>3)/3;
    uint128_t start = tmp3 + (!(tmp3&1)); // (m^4 * 2^lambda)/24
    if ((start<<1) > limit){
        throw std::runtime_error("Lamda is too large. Execution aborted. \n It should be that 2^{lambda+1} * (R.size())^4 <= 79,609,057,552,317,297,263,087,520");
    }
    std::cout << "start = " << start << std::endl;
    std::cout << "finish = " << start*2 << std::endl;
    for (uint128_t k = start; k < (start<<1); k+=2){
        if ((k - start)%5 == 0){
            std::cout << k << std::endl;
        }
        bool match = is_prime(k);
        if (match){
            uint128_t r = random_in_range_prime(k);
            return {k, r};
        }
    }
    __builtin_unreachable();
}

uint64_t find_good_prime_for_rb_tables(uint32_t N){
    uint64_t tmp = static_cast<uint64_t>(N);
    uint64_t j;
    uint64_t p = 0;
    uint64_t r = tmp;
    int k = 0;
    for (j = (tmp<<7) + 1; j <= (tmp<<8) + 1 && k < 10000; ++j){
        if (is_prime(static_cast<uint128_t>(j))){
            k++;
            if ((j%tmp) < r){
                p = j;
                r = j%tmp;
            }
        }
    }
    return p;
}

uint32_t find_prime_Bertrand(uint32_t i){ // Suppose i < 2^31
    if (i == 1){
        return 2;
    }
    uint32_t start = i + (!(i&1));
    for (uint32_t y = start; y < (i<<1); y += 2){
        if (is_prime(y)) return y;
    }
    __builtin_unreachable();
}