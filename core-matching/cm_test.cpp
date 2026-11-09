#include "core-matching.hpp"
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <string>
#include <unordered_map>
#include <chrono>

int main(){
    const uint32_t lambda_ = 2;
    const uint128_t q_prime = SUPER_PRIME;
    std::string R, T;
    std::ifstream file("../func_inversion/longtxt.txt");
    file >> R;

    R = R.substr(1764,10000);
    const uint32_t n = 1000;
    T = R.substr(30, 3*n);

    const uint32_t m = (uint32_t)R.size();
    const uint32_t query_len = 3*n;
    const uint32_t query_pos = 0;
    const uint32_t tau_ = n;

    std::cout << "m = " << m << ", |T| = " << T.size() << ", n = " << n << "\n";

    access_func access_ref  = [&R](uint32_t x) -> unsigned char { return R[x]; };
    access_func access_text = [&T](uint32_t x) -> unsigned char { return T[x]; };

    Sync_Set sync_set(lambda_, tau_);
    const uint128_t kr_base  = sync_set.base_;
    const uint128_t kr_prime = q_prime;

    const uint32_t k_sync = 9*lambda_;
    const uint32_t k = std::min(72*k_sync + 2, tau_/2);
    std::cout << "k (sparsity bound for f domain) = " << k << "\n";

    const uint32_t num_blocks = (m + n - 1)/n;
    const uint64_t N_f_64 = (uint64_t)k * num_blocks;
    if (N_f_64 > 0xFFFFFFFFull){ std::cerr << "N_f overflows uint32\n"; return 1; }
    const uint32_t N_f = (uint32_t)N_f_64;

    const uint128_t sentinel_fp = 0;

    std::function<uint128_t(uint32_t)> f_pre =
    [&, n, k, m, kr_base, kr_prime](uint32_t x) -> uint128_t {
        uint32_t b = x / k;
        uint32_t j = x - b*k;
        uint32_t block_start = b * n;
        if (block_start + 2*n > m) return 3; // before was sentinel_fp
        uint32_t win_len = std::min<uint32_t>(3*n, m - block_start);
        if (win_len < 2*n) return 2;  // before was sentinel_fp
        std::vector<uint32_t> S = sync_set.enumerate_sync(block_start, win_len, access_ref);
        uint32_t i_ref = 0;
        if (S.empty() || (S.size() <= j)){
            return mulmod((uint128_t)x, kr_base); // before was sentinel_fp
        }
        i_ref = S[j];
        if ((i_ref >= ((b*n) + n)) || ((i_ref + (n<<1)) > m)){
            return 1; // before was sentinel_fp
        }
        return phi_at_pos(i_ref, 2*n, kr_base, access_ref);
    };
    auto t0 = std::chrono::high_resolution_clock::now();
    uint128_t r = f_pre(0);
    auto t1 = std::chrono::high_resolution_clock::now();
    std::cout << "one f eval: "
            << std::chrono::duration_cast<std::chrono::microseconds>(t1-t0).count()
            << " us\n";
    std::vector<uint128_t> image_of_f;
    for (uint32_t i = 0; i < N_f; i++){
        image_of_f.push_back(f_pre(i));
    }
    std::function<uint128_t(uint32_t)> f =
    [&image_of_f](uint32_t x) -> uint128_t {
        return image_of_f[x];
    };
    const uint32_t pmax = (n + 2)/3;
    const uint32_t num_g_blocks = (m + pmax - 1)/pmax;
    const uint32_t N_g = num_g_blocks;
    const uint128_t base_powed_g = compute_pow_modulo_n_long(kr_base, (n + pmax) - 1, kr_prime);

    std::function<uint128_t(uint32_t)> g =
    [&, n, pmax, m, kr_base, kr_prime, base_powed_g](uint32_t x) -> uint128_t {
        uint32_t blk_lo = x * pmax;
        if (blk_lo + n + pmax > m) return 0;
        uint128_t roll = phi_at_pos(blk_lo, n + pmax, kr_base, access_ref);
        uint128_t best = roll;
        for (uint32_t i = blk_lo + 1; i < blk_lo + pmax; ++i){
            if (i + n + pmax > m) break;
            roll = phi_step(i, n + pmax, kr_base, base_powed_g, roll, access_ref);
            if (roll < best) best = roll;
        }
        return best;
    };

    int INRBT_f = std::max(1, (int)std::log((double)N_f));
    int INRBT_g = std::max(1, (int)std::log((double)N_g));
    std::cout << "N_f = " << N_f << " (INRBT " << INRBT_f << "), "
              << "N_g = " << N_g << " (INRBT " << INRBT_g << ")\n";
    std::cout << "building rb_main..." << std::flush;
    Mega_Rainbow_Table rb_main(INRBT_f, N_f, f);
    std::cout << " done\n";
    auto t2 = std::chrono::high_resolution_clock::now();
    uint64_t x = rb_main.invert(f(0));
    auto t3 = std::chrono::high_resolution_clock::now();
    std::cout << "one f inversion: "
            << std::chrono::duration_cast<std::chrono::microseconds>(t3-t2).count()
            << " us\n";
    uint32_t misses = 0;
    uint32_t nb_sentinels = 0;
    for (uint32_t cc = 0; cc < N_f; ++cc){
        uint128_t y = f_pre(cc);
        std::cout << "cc = " << cc << std::endl;
        if (y == 0){
            nb_sentinels++;
        }
        uint64_t x = rb_main.invert(y);
        if (x == NONE){
            misses++;
        }
    }
    std::cout << "misses = " << misses << std::endl;
    std::cout << "nb_sentinels = " << nb_sentinels << std::endl;

    std::cout << "building rb_periodic..." << std::flush;
    Mega_Rainbow_Table rb_periodic(INRBT_g, N_g, g);
    std::cout << " done\n";
    std::cout << "calling core_matching..." << std::flush;

    std::pair<uint32_t,uint32_t> res = core_matching(
        query_pos, query_len, m,
        access_text, access_ref,
        sync_set, k,
        kr_base, rb_main, rb_periodic);

    const uint32_t FAIL_I = std::numeric_limits<uint32_t>::max();
    if (res.first == FAIL_I){
        std::cout << "core_matching returned FAIL\n";
    } else {
        std::cout << "core_matching returned (i=" << res.first << ", len=" << res.second << ")\n";
    }
    return 0;
}
