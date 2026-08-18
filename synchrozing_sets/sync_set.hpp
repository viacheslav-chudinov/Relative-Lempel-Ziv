#pragma once
#include "../find_a_prime/primes.hpp"
#include "../common_utils.hpp"
#include <boost/multiprecision/cpp_int.hpp>

using access_func = std::function<char(uint32_t)>;

class Sync_Set{
private:
    const uint32_t m_;
    const uint32_t lambda_;
    const uint32_t tau_;
    const uint32_t k_;
public:
    uint128_t base_;
    uint128_t q_prime; // modulo from Karp-Rabin fingerprints. As accustomed in this project, it will be = 2^127 - 1
    bool sync(const uint32_t ind, const uint128_t kr_fingerprint, const access_func& access);
    uint32_t d_; // will be a prime
    std::vector<uint32_t> hash_coeffs;
    Sync_Set(const uint32_t m_, const uint32_t lambda_, const uint32_t tau_, const uint128_t q_prime) : m_ (m_), lambda_(lambda_),
    tau_(tau_), q_prime{q_prime}, k_(9*lambda_){
        if (tau_ < 12*k_){
            this->hash_coeffs.resize(k_, 0);
            this->d_ = 2; // Does not matter as in this case hash polynom is = 0, so d_  does not matter.
        } else{
            uint32_t start = (tau_+1)/(12*k_);
            this->d_ = find_prime_Bertrand(start);
            this->hash_coeffs.reserve(k_);
            for (uint32_t i = 0; i < k_; ++i){
                this->hash_coeffs[i] = random_in_range_prime_with_0(d_);
            }
        }
        base_ = random_in_range_prime(q_prime);
    }
    void enum_block(const uint32_t ind, const access_func& access, std::vector<uint32_t>& sync_pos);
    std::vector<uint32_t> enumerate_sync(const uint32_t ind, const uint32_t length, const access_func& access);
};

uint32_t eval_hash(uint128_t x, std::vector<uint32_t>& coeff, uint32_t d_);