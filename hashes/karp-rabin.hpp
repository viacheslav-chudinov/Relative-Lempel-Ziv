#pragma once
#include "../common_utils.hpp"

using access_func = std::function<unsigned char(uint32_t)>;

inline uint128_t phi_at_pos(const uint32_t pos, const uint32_t len, const uint128_t base, const access_func& access)
{
    uint128_t res = static_cast<uint128_t>(0);
    for (uint32_t i = 0; i < len ; ++i){
        res = addmod(mulmod(res, base), access(i+pos));
    }
    return res;
}

inline uint128_t phi_step(const uint32_t pos, const uint32_t len, uint128_t base, const uint128_t base_powed,const uint128_t previous_hash, const access_func& access)
// base_powed = base^{len - 1}, where len is obviously the length of window
{
    if (pos == 0){
        throw std::runtime_error("Bad indexing in phi_step()");
    }
    uint128_t hash = reduce((q_prime+previous_hash) - mulmod(access(pos-1),base_powed));
    return addmod(mulmod(hash, base), access(pos+len-1));
}

