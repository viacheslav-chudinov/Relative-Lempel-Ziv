#include <cmath>
#include "rainbow_table.hpp"

template <class Rainbow_Table>
class Mega_Rainbow_Table{
private:
    int NRT = 7; // An optimal number found through trial and error. Number of Rainbow Tables
    std::vector<Rainbow_Table> tables;
public:
    Mega_Rainbow_Table(uint32_t N, std::function<uint32_t(uint32_t)> func, uint64_t p){
        tables.resize(NRT);
        uint32_t cbrt_N = static_cast<uint32_t>(std::cbrt(N)); // std::cbrt(N) = N^{1/3}
        uint32_t nb_chains = 2 * cbrt_N * cbrt_N;
        uint32_t chain_length = 2 * cbrt_N; 
        for (int i = 0; i < NRT; ++i){
            tables[i] = Rainbow_Table(N, func, p, nb_chains, chain_length);
        }
    }
    uint64_t invert(uint32_t y); // Inverts a query: outputs an x such that f(x) = y
};