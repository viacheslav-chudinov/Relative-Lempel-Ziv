#pragma once
#include <algorithm>
#include <vector>
#include <functional>
#include "../common_utils.hpp"

class Rainbow_Table_3{
// Invert a function f: {0,1,..,N-1} -> {0,1,..,M-1}. AND M >= N !
// Reduction functions are : R_j : u -> ((u + c_j) mod p) mod N, where
// p is prime, c_j is a random constant in the range [1..p-1]

private:
    uint32_t N;
    std::function<uint128_t(uint32_t)> f; 
    uint128_t prime;
    uint32_t nb_chains;
    uint32_t chain_length;
    std::vector<std::pair<uint32_t, uint32_t>> table;
    std::vector<std::function<uint32_t(uint128_t)>> reductions;
    void build_reductions();
    void build_table();
public: 
    Rainbow_Table_3(uint32_t N, std::function<uint128_t(uint32_t)> func, uint128_t p, uint32_t nb_chains, 
    uint32_t chain_length) : N {N}, f {func}, prime {p}, nb_chains {nb_chains}, chain_length {chain_length}{

        this->table.resize(nb_chains, {0,0});
        this->build_reductions();
        this->build_table();
        
    }
    uint64_t invert(uint128_t y); // Outputs an x such that f(x) = y, or, if does not exist, outputs NONE = max(uint64_t)
};

