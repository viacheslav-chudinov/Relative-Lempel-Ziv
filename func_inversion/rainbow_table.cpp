#include "rainbow_table.hpp"

void Rainbow_Table_3::build_reductions(){
    // Building reduction functions for a table
    // More specifically, a chain will look like:
    // x_1 -> R_1(f(x)) =: x_2 -> R_2(f(x_2)) =: x_2 -> and so on
    reductions.clear();
    reductions.reserve(chain_length);
    for (uint j = 0; j < chain_length; j++){
        uint128_t c_j = random_in_range_prime(this->prime);
        reductions.push_back([c_j, this](uint128_t h){
            uint128_t mod_p = (h + c_j)%this->prime;
            return static_cast<uint32_t>(mod_p % N);
        });
    }
}

void Rainbow_Table_3::build_table(){
    uint32_t k,l;
    for (k = 0; k < nb_chains; k++){
        uint32_t x0 = random_in_range(N); 
        uint32_t x = x0;
        for (l = 0; l < chain_length; ++l){
            x = reductions[l](f(x));
        }
        table[k] = {x0, x};
    }

    std::sort(table.begin(), table.end(), 
    [](const std::pair<uint32_t,uint32_t>& a, const std::pair<uint32_t,uint32_t>& b){ return a.second < b.second; });
}


uint64_t Rainbow_Table_3::invert(uint128_t y){
    for (long int c = chain_length - 1; c >= 0; --c){
        uint32_t endpoint = reductions[c](y);
        for (long int k = c + 1; k < (long int)chain_length; ++k){
            endpoint = reductions[k](f(endpoint));
        }
        auto it = std::lower_bound(table.begin(), table.end(), endpoint,
            [](const std::pair<uint32_t,uint32_t>& a, uint32_t val){ return a.second < val; });
        while (it != table.end() && it->second == endpoint){
            uint32_t x = it->first;
            for (long int k = 0; k < (long int)chain_length; ++k){
                if (f(x) == y){
                    return x;       
                }
                x = reductions[k](f(x));
            }
            ++it;                   
        }
    }
    return NONE;
}