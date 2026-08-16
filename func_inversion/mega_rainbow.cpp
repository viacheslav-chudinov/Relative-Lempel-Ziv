#include "mega_rainbow.hpp"

uint64_t Mega_Rainbow_Table::invert(unsigned __int128 y){
    for (uint32_t i = 0; i < tables.size(); i++){
        uint64_t x = tables[i]->invert(y);
        if (x != NONE){
            return x;
        }
    }
    return NONE;
}

bool Mega_Rainbow_Table::check_completeness(){
    for (uint32_t k = 0; k < N; k++){
        unsigned __int128 y = f(k);
        uint64_t x = this->invert(y);
        if (x == NONE){
            return false;
        }
    }
    return true;
}

bool Mega_Rainbow_Table::try_complete(){
    uint32_t old_size = tables.size();
    uint32_t bound = std::log(N);
    while ((tables.size() - old_size) < bound){
        Rainbow_Table_3* rb = new Rainbow_Table_3(N, f, prime, nb_chains, chain_length);
        tables.push_back(rb);
        if (this->check_completeness()){
            return true;
        }
    }
    return false;
}

void Mega_Rainbow_Table::pretty_add_one_table(){
    Rainbow_Table_3* rb = new Rainbow_Table_3(N, f, prime, nb_chains, chain_length);
    tables.push_back(rb);
}