#include <cmath>
#include <set>
#include <unordered_set>
#include "rainbow_table.hpp"

// Used specifically for Karp-Rabin fingerprints
class Mega_Rainbow_Table{
private:
    int INRBT; // Initial Number of RainBow Tables
    uint32_t N;
    // std::vector<Rainbow_Table_3*> tables;
    uint128_t prime;
    uint32_t nb_chains;
    uint32_t chain_length;
public:
    std::function<uint128_t(uint32_t)> f;
    std::vector<Rainbow_Table_3*> tables;
    Mega_Rainbow_Table(int INRBT, uint32_t N, std::function<uint128_t(uint32_t)> func, uint128_t p) : INRBT {INRBT}, N {N}, 
    f{func}{
        tables.resize(INRBT);
        uint32_t cbrt_N = static_cast<uint32_t>(std::cbrt(N)); // std::cbrt(N) = N^{1/3}
        this->nb_chains = 2 * cbrt_N * cbrt_N;
        this->chain_length = 2 * cbrt_N; 
        for (int i = 0; i < INRBT; ++i){
            tables[i] = new Rainbow_Table_3(N, func, p, this->nb_chains, this->chain_length);
        }
    }
    ~Mega_Rainbow_Table(){
        for (auto* t : tables) delete t;
    }
    bool check_completeness();
    bool try_complete();
    void pretty_add_one_table();
    uint64_t invert(uint128_t y); // Inverts a query: outputs an x such that f(x) = y
};