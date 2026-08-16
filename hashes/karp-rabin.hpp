#include "../common_utils.hpp"

using access_func = std::function<unsigned char(uint32_t)>;
using boost::multiprecision::uint256_t;

uint128_t phi_at_pos(const uint32_t pos, const uint32_t len, const uint128_t base, const uint128_t prime, const access_func& access)
{
    uint256_t res = static_cast<uint256_t>(0);
    uint256_t P = static_cast<uint256_t>(prime);
    for (uint32_t i = 0; i < len ; ++i){
        res = (((res * base)%P) + access(i + pos))%P;
    }
    return static_cast<uint128_t>(res);
}

uint128_t phi_step(const uint32_t pos, const uint32_t len, uint128_t base, const uint128_t base_powed,const uint128_t previous_hash, const uint128_t prime, const access_func& access)
// base_powed = base^{len - 1}, where len is obviously the length of window
{
    if (pos == 0){
        throw std::runtime_error("Bad indexing in phi_step()");
    }
    uint256_t hash = static_cast<uint256_t>(previous_hash);
    uint256_t P = static_cast<uint256_t>(prime);
    uint256_t base_powed_256 = static_cast<uint256_t>(base_powed);
    hash = ((P+hash) - ((access(pos-1)*base_powed_256)%P))%P;
    hash = ((hash*(uint256_t)base)%P + (access(pos+len-1)))%P;
    return static_cast<uint128_t>(hash);
}

