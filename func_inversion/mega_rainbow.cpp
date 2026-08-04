#include "mega_rainbow.hpp"

template <class Rainbow_Table>
uint64_t Mega_Rainbow_Table<Rainbow_Table>::invert(uint32_t y){
    for (int i = 0; i < NRT; i++){
        uint64_t x = tables[i].invert(y);
        if (x != NONE){
            return x;
        }
    }
    return NONE;
}