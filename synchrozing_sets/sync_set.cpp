#include "sync_set.hpp"
#include "../compute_period.hpp"

uint32_t eval_hash(uint128_t x, std::vector<uint32_t>& coeffs, uint32_t d_){
    using boost::multiprecision::uint256_t;
    uint256_t res = static_cast<uint256_t>(0);
    uint256_t P = static_cast<uint256_t>(d_);
    for (long int i = coeffs.size() - 1; i >= 0; --i)
        res = (((res * x)%P) + coeffs[i])%P;
    return static_cast<uint128_t>(res);
}

bool Sync_Set::sync(const uint32_t ind, const uint128_t kr_fingerprint, const access_func& access){
    uint32_t per_plus = compute_period(ind + 1, tau_<<1, access);
    uint32_t per_minus = compute_period(ind, (tau_<<1) - 1, access);
    if (per_plus > tau_/3 && per_minus > tau_/3){
        return eval_hash(kr_fingerprint, hash_coeffs, d_) == 0;
    }
    uint32_t per = compute_period(ind, tau_<<1, access);
    if (per > tau_/3){
        return true;
    } 
    return false;
}