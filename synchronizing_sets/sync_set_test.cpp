#include "sync_set.hpp"
#include <iostream>
#include <fstream>
#include "../hashes/karp-rabin.hpp"


int main(){
    const uint32_t lambda_ = 6; 
    const uint32_t tau_ = 500000;
    Sync_Set sync_set(lambda_, tau_);
    std::string R;
    std::ifstream file("../func_inversion/longtxt.txt");
    std::string text;
    file >> R;
    std::cout << "R's size is " << R.size() << std::endl;
    access_func access = [&R](uint32_t x) -> unsigned char{
        return R[x];
    };
    const uint32_t ind = 0;
    const uint32_t length = R.size();
    std::vector<uint32_t> sync_pos = sync_set.enumerate_sync(ind, length, access);
    std::cout << "sync pos.size() = " << sync_pos.size() << std::endl;
    // for (uint32_t i = 0; i < sync_pos.size(); ++i){
    //     std::cout << "i-th, i = " << i + 1 << " synchronizing position is " << sync_pos[i] << std::endl;
    // }
    std::cout << "d = " << sync_set.d_ << "\n";
    // std::cout << "coeffs: ";
    // for (auto c : sync_set.hash_coeffs) std::cout << c << " ";
    // std::cout << "\n";
    // std::cout << "base = " << (uint64_t)sync_set.base_ << "\n";  // low bits enough to compare
    std::cout << "Sparsity and density status is " << sync_set.is_sparse_and_dense(ind, length, access) << std::endl;
    // int hits = 0;
    // for (uint32_t p = 0; p + 2*tau_ <= length; p += 500) {
    //     uint128_t fp = phi_at_pos(p, 2*tau_, sync_set.base_, q_prime, access);
    //     if (eval_hash(fp, sync_set.hash_coeffs, sync_set.q_prime, sync_set.d_) == 0) hits++;
    // }
    // std::cout << "raw hash hits: " << hits << "\n";
    return 0;
}
