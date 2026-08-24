#include "../synchronizing_sets/sync_set.hpp"
#include "../compute_period.hpp"
#include "../hashes/karp-rabin.hpp"
#include "../common_utils.hpp"
#include "../func_inversion/mega_rainbow.hpp"
#include <vector>

std::pair<uint32_t, uint32_t> core_matching(
    uint32_t query_pos, uint32_t query_len, uint32_t ref_length,
    access_func& access_text, access_func& access_ref,
    Sync_Set& sync_set, uint32_t k,
    uint128_t kr_base, Mega_Rainbow_Table& rb_table_main,
    Mega_Rainbow_Table& rb_table_periodic_case);