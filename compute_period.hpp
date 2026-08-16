#include "common_utils.hpp"
#include <boost/multiprecision/cpp_int.hpp>

using access_func = std::function<char(uint32_t)>;


uint32_t compute_period(const uint32_t pos, const uint32_t len, const access_func& access){
    if (len <= 1){
        return len; 
    }
    uint32_t d = (len + 1) / 2;
    auto w_at = [pos, &access](uint32_t k) -> unsigned char {
        return access(pos + k);
    };

    auto UpdateMS = [&w_at](uint32_t l, uint32_t s, uint32_t p) -> std::tuple<uint32_t,uint32_t,uint32_t> {
        if (l == 0){
            return {1, 0, 1};
        }
        uint32_t i = l;
        while (i < l + 1){
            if (w_at(i - p) < w_at(i)){
                i = i - (i >= s ? (i - s) % p : p - ((s - i) % p));
                s = i;
                p = 1;
            } else if (w_at(i - p) > w_at(i)){
                p = i - s + 1;
            }
            i++;
        }
        return {l + 1, s, p};
    };

    uint32_t i = 1, p = 0, s = 0, l = 0;
    const uint32_t m = d;      
    const uint32_t n = len;
    while (i + d <= n){    
        while (i + l < n && l < m && w_at(i + l) == w_at(l)){
            std::tie(l, s, p) = UpdateMS(l, s, p);
        }
        if (l == m){
            uint32_t period_candidate = i;
            bool is_period = true;
            for (uint32_t j = 0; j + period_candidate < len; ++j){
                if (access(pos + j) != access(pos + j + period_candidate)){
                    is_period = false;
                    break;
                }
            }
            return is_period ? period_candidate : len;
        }

        bool s_match = true;
        for (uint32_t b = 0; b < s; ++b){
            if (w_at(b) != w_at(p + b)){
                s_match = false;
                break;
            }
        }
        if (p > 0 && p <= l / 3 && s_match){
            i += p;
            l -= p;
        } else {
            i += l / 3 + 1;
            l = 0; s = 0; p = 0;
        }
    }
    return len; 
}