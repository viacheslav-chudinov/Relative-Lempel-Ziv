#include "../synchrozing_sets/sync_set.hpp"
#include "../compute_period.hpp"
#include "../hashes/karp-rabin.hpp"
#include "../common_utils.hpp"
#include "../func_inversion/mega_rainbow.hpp"
#include <vector>

std::pair<uint32_t, uint32_t> core_matching(
    uint32_t query_pos, uint32_t query_len, uint32_t ref_length,
    access_func& access_text, access_func& access_ref,
    Sync_Set& sync_set, uint32_t k,
    uint128_t kr_base, uint128_t kr_prime,
    Mega_Rainbow_Table& rb_table_main,
    Mega_Rainbow_Table& rb_table_periodic_case)
{
    const uint32_t FAIL_I = std::numeric_limits<uint32_t>::max();
    const auto FAIL = std::make_pair(FAIL_I, (uint32_t)0);
    const uint32_t n = query_len / 3;

    std::vector<uint32_t> Sq = sync_set.enumerate_sync(query_pos, query_len, access_text);
    bool found = !Sq.empty();
    uint32_t jj_found = 0;
    uint128_t y = 0;
    if (found){
        jj_found = Sq.front() - query_pos;
        y = phi_at_pos(query_pos + jj_found, (n<<1) , kr_base, kr_prime, access_text);
        uint64_t x = rb_table_main.invert(y);
        if (x == NONE){
            return FAIL;
        }
        uint32_t b = (uint32_t)(x / k);               
        uint32_t j = (uint32_t)(x - (uint64_t)b * k);
        uint32_t candidate_start = b * n;
        if (candidate_start + (n<<1) > ref_length){ 
            throw std::runtime_error("Candidate start is off");
        }
        uint32_t candidate_len = std::min<uint32_t>(query_len, ref_length - candidate_start);
        std::vector<uint32_t> S = sync_set.enumerate_sync(candidate_start, candidate_len, access_ref);

        uint32_t i_ref = FAIL_I;
        if (S.empty() || (S.size() <= j)){
            return FAIL;
        }
        i_ref = S[j];
        if ((i_ref >= (candidate_start + n)) || ((i_ref + (n<<1)) > ref_length)){
            return FAIL;
        }
        for (uint32_t t = 0; t < (n<<1); ++t){
            if (access_ref(i_ref + t) != access_text(query_pos + jj_found + t)){
                return FAIL;
            }
        }
        uint32_t occ_start = i_ref + n - jj_found;
        return {occ_start, n};
    }
    const uint32_t pmax = (n + 2) / 3;
    uint32_t j_prime = n - pmax;
    uint128_t base_powed_per = compute_pow_modulo_n_long(kr_base, (n + pmax) - 1, kr_prime);
    uint128_t roll_kr = phi_at_pos(query_pos + (n - pmax), n + pmax, kr_base, kr_prime, access_text);
    uint128_t y_per = roll_kr; 
    for (uint32_t j = n - pmax + 1; j <= n - 1; ++j){

        roll_kr = phi_step(query_pos + j, n + pmax, kr_base, base_powed_per,
                roll_kr, kr_prime, access_text);

        if (roll_kr < y_per){ 
            y_per = roll_kr; 
            j_prime = j; 
        }
    }
    uint64_t x = rb_table_periodic_case.invert(y_per);
    if (x == NONE){
        return FAIL;
    }
    uint32_t left_start = x * pmax; // checking the interval [x*p_max...x*p_max+p_max)
    uint32_t i_ref = FAIL_I; 
    uint128_t best_kr = 0;
    if (left_start + n + pmax <= ref_length){
        uint128_t roll_r = phi_at_pos(left_start, n + pmax, kr_base, kr_prime, access_ref);
        best_kr = roll_r; 
        i_ref = left_start;
        for (uint32_t i = left_start + 1; i < left_start + pmax; ++i) {
            if (i + n + pmax > ref_length) break;

            roll_r = phi_step(i, n + pmax, kr_base, base_powed_per,
                              roll_r, kr_prime, access_ref);

            if (roll_r < best_kr){
                best_kr = roll_r; 
                i_ref = i; 
            }
        }
    }
    if (i_ref == FAIL_I){
        return FAIL;
    }
    for (uint32_t t = 0; t < n + pmax; ++t){
        if (access_ref(i_ref + t) != access_text(query_pos + j_prime + t)){
            return FAIL;
        }
    }
    uint32_t occ_start = i_ref + n - j_prime;
    return {occ_start, n};
}