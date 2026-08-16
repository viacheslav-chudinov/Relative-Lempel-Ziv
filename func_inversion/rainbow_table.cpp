#include "rainbow_table.hpp"


void Rainbow_Table_1::build_reductions(){
    // Building reduction functions for a table
    // More specifically, a chain will look like:
    // x_1 -> R_1(f(x)) =: x_2 -> R_2(f(x_2)) =: x_2 -> and so on
    reductions.clear();
    reductions.reserve(chain_length);
    static std::random_device rd;
    static std::mt19937_64 gen(rd());
    std::uniform_int_distribution<uint64_t> dist(1, prime - 1);
    for (uint j = 0; j < chain_length; j++){
        uint64_t c_j = random_in_range_prime(prime);
        reductions.push_back([c_j, this](uint64_t h) {
            uint128_t tmp = (uint128_t)(h + c_j);
            uint64_t mod_p = (uint64_t)(tmp % prime);
            return (uint32_t)(mod_p % N);
        });
    }
}

void Rainbow_Table_2::build_reductions(){
    // Building reduction functions for a table
    // More specifically, a chain will look like:
    // x_1 -> R_1(f(x)) =: x_2 -> R_2(f(x_2)) =: x_2 -> and so on
    reductions.clear();
    reductions.reserve(chain_length);
    static std::random_device rd;
    static std::mt19937_64 gen(rd());
    std::uniform_int_distribution<uint64_t> dist(1, prime - 1);
    for (uint j = 0; j < chain_length; j++){
        uint64_t a_j = random_in_range_prime(prime);
        uint64_t b_j = random_in_range_prime(prime);
        reductions.push_back([a_j, b_j, this](uint32_t h) {
            uint128_t tmp = (uint128_t)a_j * h + b_j;
            uint64_t mod_p = (uint64_t)(tmp % prime);
            return (uint32_t)(mod_p % N);
        });
    }
}

void Rainbow_Table_3::build_reductions(){
    // Building reduction functions for a table
    // More specifically, a chain will look like:
    // x_1 -> R_1(f(x)) =: x_2 -> R_2(f(x_2)) =: x_2 -> and so on
    reductions.clear();
    reductions.reserve(chain_length);
    static std::random_device rd;
    static std::mt19937_64 gen(rd());
    std::uniform_int_distribution<uint64_t> dist(1, prime - 1);
    for (uint j = 0; j < chain_length; j++){
        uint64_t c_j = random_in_range_prime(prime);
        reductions.push_back([c_j, this](uint64_t h) {
            uint128_t tmp = (uint128_t)(h + c_j);
            uint64_t mod_p = (uint64_t)(tmp % prime);
            return (uint32_t)(mod_p % N);
        });
    }
}

void Rainbow_Table_0::build_table(){
    uint32_t k,l;
    for (k = 0; k < nb_chains; k++){
        uint32_t x0 = random_in_range(N); 
        uint32_t x = x0;
        for (l = 0; l < chain_length; ++l){
            x = f(x);
        }
        table[k] = {x0, x};
    }
    
    std::sort(table.begin(), table.end(), 
    [](const std::pair<uint32_t,uint32_t>& a, const std::pair<uint32_t,uint32_t>& b){ return a.second < b.second; });
}

void Rainbow_Table_1::build_table(){
    uint32_t k,l;
    for (k = 0; k < nb_chains; k++){
        uint32_t x0 = random_in_range(N); 
        uint32_t x = x0;
        for (l = 0; l < chain_length; ++l){
            x = f(reductions[l](x));
        }
        table[k] = {x0, x};
    }

    std::sort(table.begin(), table.end(), 
    [](const std::pair<uint32_t,uint32_t>& a, const std::pair<uint32_t,uint32_t>& b){ return a.second < b.second; });
}

void Rainbow_Table_2::build_table(){
    uint32_t k,l;
    for (k = 0; k < nb_chains; k++){
        uint32_t x0 = random_in_range(N); 
        uint32_t x = x0;
        for (l = 0; l < chain_length; ++l){
            x = f(reductions[l](x));
        }
        table[k] = {x0, x};
    }

    std::sort(table.begin(), table.end(), 
    [](const std::pair<uint32_t,uint32_t>& a, const std::pair<uint32_t,uint32_t>& b){ return a.second < b.second; });
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


uint64_t Rainbow_Table_0::invert(uint32_t y){
    for (long int j = chain_length; j >= 0 ; --j){
        uint32_t multi_y = y;
        for (long int k = j; k < chain_length; ++k){
            multi_y = f(multi_y);
        }
        auto it = std::lower_bound(table.begin(), table.end(), 
        multi_y,[](const std::pair<uint32_t,uint32_t> &a, uint32_t val)
        {
            return a.second < val;
        });
        while (it != table.end() && it->second == multi_y){
            uint32_t x = it->first;
            for (long int k = 0; k < chain_length; ++k){
                if ( f(x) == y ){
                    return x;
                }
                x = f(x);
            }
            ++it; 
        }
    }
    return NONE;
}
uint64_t Rainbow_Table_1::invert(uint32_t y){
    for (long int j = chain_length; j >= 0 ; --j){
        uint32_t multi_y = y;
        for (long int k = j; k < chain_length; ++k){
            multi_y = f(reductions[k](multi_y));
        }
        auto it = std::lower_bound(table.begin(), table.end(), 
        multi_y,[](const std::pair<uint32_t,uint32_t> &a, uint32_t val)
        {
            return a.second < val;
        });
        while (it != table.end() && it->second == multi_y){
            uint32_t x = it->first;
            uint32_t tmp;
            for (long int k = 0; k < chain_length; ++k){
                tmp = f(reductions[k](x));
                if ( tmp == y ){
                    return reductions[k](x);
                }
                x = f(reductions[k](x));
            }
            ++it; 
        }
    }
    return NONE;
}

uint64_t Rainbow_Table_2::invert(uint32_t y){
    for (long int j = chain_length; j >= 0 ; --j){
        uint32_t multi_y = y;
        for (long int k = j; k < chain_length; ++k){
            multi_y = f(reductions[k](multi_y));
        }
        auto it = std::lower_bound(table.begin(), table.end(), 
        multi_y,[](const std::pair<uint32_t,uint32_t> &a, uint32_t val)
        {
            return a.second < val;
        });
        while (it != table.end() && it->second == multi_y){
            uint32_t x = it->first;
            uint32_t tmp;
            for (long int k = 0; k < chain_length; ++k){
                tmp = reductions[k](x);
                if ( f(tmp) == y ){
                    return tmp;
                }
                x = f(tmp);
            }
            ++it; 
        }
    }
    return NONE;
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