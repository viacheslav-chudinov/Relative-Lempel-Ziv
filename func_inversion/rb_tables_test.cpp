#include <iostream>
#include <fstream>
#include <set>
#include <unordered_set>
#include "rainbow_table.hpp"
#include "../find_a_prime/primes.hpp"
#include <boost/multiprecision/cpp_int.hpp>
#include "mega_rainbow.hpp"

uint32_t func_1(uint32_t x){ // x <= 50
    if (x%2 == 0){
        return (x+2)%52; // Encodes cycle 0->2->...->50->0
    } else if(x%4 == 1 && x <= 33){
        return (x+4)%36; // Encodes cycle 1->5->...->33->1
    } else if (x%4 == 3 && x <= 39){
        return (x+4)%40; // Encodes cycle 3->7->...->39->3
    } else if (x == 37){
        return 43;
    } else if (x == 43){
        return 49;
    } else if (x == 49){
        return 41;
    } else if (x == 41){
        return 37;
    } else if (x == 45){
        return 45;
    } else{ // x = 47
        return 47;
    }
}

uint32_t compute_polynom_modulo_n_unsigned(uint32_t x, std::vector<std::pair<uint32_t, uint32_t>>& coeffs, uint32_t n){
    uint32_t k = 0;
    uint64_t ans = 0;
    while(k < coeffs.size()){
        auto [pow, coef] = coeffs[k];
        ans = (ans + (coef*compute_pow_modulo_n(x, pow, n))%n)%n;
        k++;
    }
    return ans;
}


uint32_t func_3(uint32_t x){ // x <= 10_000
    static const uint32_t p = 42043; // A prime number
    static std::random_device rd;
    static std::mt19937_64 gen(rd());
    static std::uniform_int_distribution<uint32_t> dist(0, p-1);
    static const int deg = 250;
    static std::vector<std::pair<uint32_t, uint32_t>> coeffs = [](){
        std::random_device rd;
        std::mt19937_64 gen(rd());
        std::uniform_int_distribution<uint32_t> dist(0, p - 1);
        std::vector<std::pair<uint32_t, uint32_t>> c;
        for (int r = deg; r >= 0; r--){
            c.push_back({r, dist(gen)});
        }
        return c;
    }();
    return compute_polynom_modulo_n_unsigned(x, coeffs, p);
}

uint128_t evalute_poly(uint128_t x, const std::vector<uint32_t>& coeffs, uint128_t p){
    using boost::multiprecision::uint256_t;
    uint256_t res = static_cast<uint256_t>(0);
    uint256_t P = static_cast<uint256_t>(p);
    for (long int i = coeffs.size() - 1; i >= 0; --i)
        res = (((res * x)%P) + coeffs[i])%P;
    return static_cast<uint128_t>(res);
}

uint32_t func_4(uint32_t x){ // x <= 42043
    static const uint32_t p = 42043; // A prime number
    static std::random_device rd;
    static std::mt19937_64 gen(rd());
    static std::uniform_int_distribution<uint32_t> dist(0, p-1);
    static const int deg = 1000;
    static std::vector<uint32_t> coeffs = [](){
        std::random_device rd;
        std::mt19937_64 gen(rd());
        std::uniform_int_distribution<uint32_t> dist(0, p - 1);
        std::vector<uint32_t> c;
        for (int r = deg; r >= 0; r--){
            c.push_back(dist(gen));
        }
        return c;
    }();
    return evalute_poly(x, coeffs, p) + 1;
}

std::vector<uint32_t> func_5_KR(){ 
    uint32_t p = 10781; // A prime number
    uint32_t base = random_in_range_prime(p);
    std::ifstream file("covid_rna_1000.txt");
    std::string rna;
    file >> rna;
    uint32_t window = 40;
    std::vector<uint32_t> vals;
    std::vector<uint32_t> coeffs;
    for (int x = window - 1; x >= 0; --x){
        coeffs.push_back((uint32_t)rna[x]);
    }
    uint32_t hash = evalute_poly(base, coeffs, p);
    uint32_t base_pow_window = compute_pow_modulo_n(base, window - 1, p);
    vals.push_back(hash);
    for (uint32_t x = window; x < rna.size(); ++x){
        hash = (p + (hash - ((uint32_t)(unsigned char)rna[x- window])*base_pow_window)%p)%p;
        hash = ((hash*base)%p + ((uint32_t)(unsigned char)rna[x]))%p;
        vals.push_back(hash);
    }
    return vals;
}

uint32_t func_5(uint32_t x){ // x <= rna.size() - window
    static const std::vector<uint32_t> vals = func_5_KR();
    return vals[x];
}

uint32_t func_6_prep(){
    std::ifstream file("longtxt.txt");
    std::string text;
    file >> text;
    return text.size();
}

std::vector<uint128_t> func_6_KR(){
    using boost::multiprecision::uint256_t;
    uint128_t p = SUPER_PRIME;
    uint128_t base = random_in_range_prime(SUPER_PRIME);
    std::ifstream file("longtxt.txt");
    std::string text;
    file >> text;
    uint32_t window = 40001;
    std::vector<uint128_t> vals;
    std::vector<uint32_t> coeffs(window, 0);
    for (int x = window - 1; x >= 0; --x){
        coeffs.push_back((uint32_t)text[x]);
    }
    uint256_t hash = uint256_t(evalute_poly(base, coeffs, p));
    uint256_t base_pow_window = uint256_t(compute_pow_modulo_n_long(base, window - 1, p));
    vals.push_back(static_cast<uint128_t>(hash));
    for (uint32_t x = window; x < text.size(); ++x){
        hash = ((p+hash) - (((unsigned char)text[x- window])*base_pow_window)%p)%p;
        hash = ((hash*base)%p + ((uint128_t)(unsigned char)text[x]))%p;
        vals.push_back(static_cast<uint128_t>(hash));
    }
    return vals;
}

uint128_t func_6(uint32_t x){
    static const std::vector<uint128_t> vals = func_6_KR();
    return vals[x];
}

int main(){
    // Possible tests:
    // UNCOMMENT ONLY 1 FUNCTION TEST (YES, I WAS LAZE TO RENAME VARIABLES)

    // std::cout << " ---- Function 1 ---- " << std::endl;
    // Rainbow_Table_0 rb00(51, func_1, 16,4);
    // Rainbow_Table_0 rb01(51, func_1, 16,4);
    // Rainbow_Table_1 rb10(51, func_1, 131, 16,4);
    // Rainbow_Table_1 rb11(51, func_1, 139, 16,4);
    // Rainbow_Table_2 rb20(51, func_1, 257, 16, 4);
    // Rainbow_Table_2 rb21(51, func_1, 257, 16, 4);
    // int k;
    // int rb00_misses = 0; 
    // int rb01_misses = 0; 
    // int rb10_misses = 0; 
    // int rb11_misses = 0; 
    // int rb20_misses = 0; 
    // int rb21_misses = 0; 
    // uint64_t x;
    // for (k = 0; k <= 50; k++){
    //     std::cout << "Inverting with rb00" << std::endl;
    //     x = rb00.invert(k);
    //     if (x == NONE){
    //         rb00_misses++;
    //         std::cout << "Ups, could not invert on rb00" << std::endl;
    //     }
    //     std::cout << "Inverting with rb01" << std::endl;
    //     x = rb01.invert(k);
    //     if (x == NONE){
    //         rb01_misses++;
    //         std::cout << "Ups, could not invert on rb01" << std::endl;
    //     }
    //     std::cout << "Inverting with rb10" << std::endl;
    //     x = rb10.invert(k);
    //     if (x == NONE){
    //         rb10_misses++;
    //         std::cout << "Ups, could not invert on rb10" << std::endl;
    //     }
    //     std::cout << "Inverting with rb11" << std::endl;
    //     x = rb11.invert(k);
    //     if (x == NONE){
    //         rb11_misses++;
    //         std::cout << "Ups, could not invert on rb11" << std::endl;
    //     }
    //     std::cout << "Inverting with rb20" << std::endl;
    //     x = rb20.invert(k);
    //     if (x == NONE){
    //         rb20_misses++;
    //         std::cout << "Ups, could not invert on rb20" << std::endl;
    //     }
    //     std::cout << "Inverting with rb21" << std::endl;
    //     x = rb21.invert(k);
    //     if (x == NONE){
    //         rb21_misses++;
    //         std::cout << "Ups, could not invert on rb21" << std::endl;
    //     }
    // }
    // std::cout << "00-th table got " << rb00_misses << " misses" << std::endl;
    // std::cout << "01-th table got " << rb01_misses << " misses" << std::endl;
    // std::cout << "10-st table got " << rb10_misses << " misses" << std::endl;
    // std::cout << "11-th table got " << rb11_misses << " misses" << std::endl;
    // std::cout << "20-nd table got " << rb20_misses << " misses" << std::endl;
    // std::cout << "21-th table got " << rb21_misses << " misses" << std::endl;
    // std::cout << std::endl;
    // std::cout << std::endl;




    // std::cout << " ---- Function 2 ---- " << std::endl;

    // Rainbow_Table_0 rb00(5903, func_2, 600, 16);
    // Rainbow_Table_0 rb01(5903, func_2, 300, 32);
    // Rainbow_Table_1 rb10(5903, func_2, 9431, 600, 10);
    // Rainbow_Table_1 rb11(5903, func_2, 9431, 400, 30);
    // Rainbow_Table_2 rb20(5903, func_2, 9431, 600, 10);
    // Rainbow_Table_2 rb21(5903, func_2, 9431, 500, 50);
    // std::set<uint32_t> image;
    // int k;
    // for (k = 1; k <= 5000; ++k){
    //     image.insert(func_2(k));
    // }
    // std::cout << "We have " << image.size() << " elements in the image" << std::endl;
    // int rb00_misses = 0; 
    // int rb01_misses = 0; 
    // int rb10_misses = 0; 
    // int rb11_misses = 0; 
    // int rb20_misses = 0; 
    // int rb21_misses = 0; 
    // uint64_t x;
    // int u = 0;
    // for (uint32_t y : image){
    //     if(u%20 == 0){
    //         std::cout << "u = " << u << std::endl;
    //     }
    //     u++;
    //     // std::cout << "Inverting with rb00" << std::endl;
    //     x = rb00.invert(y);
    //     if (x == NONE){
    //         rb00_misses++;
    //         // std::cout << "Ups, could not invert on rb00" << std::endl;
    //     }
    //     // std::cout << "Inverting with rb01" << std::endl;
    //     x = rb01.invert(y);
    //     if (x == NONE){
    //         rb01_misses++;
    //         // std::cout << "Ups, could not invert on rb01" << std::endl;
    //     }
    //     // std::cout << "Inverting with rb10" << std::endl;
    //     x = rb10.invert(y);
    //     if (x == NONE){
    //         rb10_misses++;
    //         // std::cout << "Ups, could not invert on rb10" << std::endl;
    //     }
    //     // std::cout << "Inverting with rb11" << std::endl;
    //     x = rb11.invert(y);
    //     if (x == NONE){
    //         rb11_misses++;
    //         // std::cout << "Ups, could not invert on rb11" << std::endl;
    //     }
    //     // std::cout << "Inverting with rb20" << std::endl;
    //     x = rb20.invert(y);
    //     if (x == NONE){
    //         rb20_misses++;
    //         // std::cout << "Ups, could not invert on rb20" << std::endl;
    //     }
    //     // std::cout << "Inverting with rb21" << std::endl;
    //     x = rb21.invert(y);
    //     if (x == NONE){
    //         rb21_misses++;
    //         // std::cout << "Ups, could not invert on rb21" << std::endl;
    //     }
    // }
    // std::cout << "00-th table got " << rb00_misses << " misses" << std::endl;
    // std::cout << "01-th table got " << rb01_misses << " misses" << std::endl;
    // std::cout << "10-st table got " << rb10_misses << " misses" << std::endl;
    // std::cout << "11-th table got " << rb11_misses << " misses" << std::endl;
    // std::cout << "20-nd table got " << rb20_misses << " misses" << std::endl;
    // std::cout << "21-th table got " << rb21_misses << " misses" << std::endl;
    // std::cout << std::endl;
    // std::cout << std::endl;




    // std::cout << " ---- Function 3 ---- " << std::endl;
    // uint32_t p = 42043;
    // Rainbow_Table_1 rb10(p, func_3, 87917, 2000, 100);
    // Rainbow_Table_1 rb11(p, func_3, 85853, 2000, 40);
    // Rainbow_Table_2 rb20(p, func_3, 87917, 2000, 100);
    // Rainbow_Table_2 rb21(p, func_3, 97729, 2000, 40);
    // std::unordered_set<uint32_t> image;
    // uint32_t k;
    // for (k = 1; k <= p; ++k){
    //     image.insert(func_3(k));
    // }
    // std::cout << "We have " << image.size() << " elements in the image" << std::endl;
    // int rb10_misses = 0; 
    // int rb11_misses = 0; 
    // int rb20_misses = 0; 
    // int rb21_misses = 0; 
    // uint64_t x;
    // int u = 0;
    // for (uint32_t y : image){
    //     u++;
    //     if (u == 400){
    //         std::cout << "u = " << u << std::endl;
    //         std::cout << "y = " << y << std::endl;
    //         std::cout << "10-st table got " << rb10_misses << " misses" << std::endl;
    //         std::cout << "11-th table got " << rb11_misses << " misses" << std::endl;
    //         std::cout << "20-nd table got " << rb20_misses << " misses" << std::endl;
    //         std::cout << "21-th table got " << rb21_misses << " misses" << std::endl;
    //         std::cout << std::endl;
    //         std::cout << std::endl;
    //         break;
    //     }
    //     std::cout << "u = " << u << std::endl;
    //     std::cout << "y = " << y << std::endl;
    //     x = rb10.invert(y);
    //     if (x == NONE){
    //         rb10_misses++;
    //     }
    //     x = rb11.invert(y);
    //     if (x == NONE){
    //         rb11_misses++;
    //     }
    //     x = rb20.invert(y);
    //     if (x == NONE){
    //         rb20_misses++;
    //     }
    //     x = rb21.invert(y);
    //     if (x == NONE){
    //         rb21_misses++;
    //     }
    // }
    // std::cout << "10-st table got " << rb10_misses << " misses" << std::endl;
    // std::cout << "11-th table got " << rb11_misses << " misses" << std::endl;
    // std::cout << "20-nd table got " << rb20_misses << " misses" << std::endl;
    // std::cout << "21-th table got " << rb21_misses << " misses" << std::endl;
    // std::cout << std::endl;
    // std::cout << std::endl;
    // std::cout << " ---- Function 4 ---- " << std::endl;
    // uint32_t p = 42043;
    // // Rainbow_Table_1 rb10(p, func_4, 87917, 600, 600);
    // Rainbow_Table_1 rb11(p, func_4, 87917, 1000, 400);
    // std::unordered_set<uint32_t> image;
    // uint32_t k;
    // for (k = 1; k <= p; ++k){
    //     image.insert(func_4(k));
    // }
    // std::cout << "We have " << image.size() << " elements in the image" << std::endl;
    // // int rb10_misses = 0; 
    // int rb11_misses = 0; 
    // uint64_t x;
    // int u = 0;
    // for (uint32_t y : image){
    //     u++;
    //     if (u == 100){
    //         std::cout << "u = " << u << std::endl;
    //         std::cout << "y = " << y << std::endl;
    //         // std::cout << "10 table got " << rb10_misses << " misses" << std::endl;
    //         std::cout << "11 table got " << rb11_misses << " misses" << std::endl;
    //         std::cout << std::endl;
    //         std::cout << std::endl;
    //         break;
    //     }
    //     std::cout << "u = " << u << std::endl;
    //     // x = rb10.invert(y);
    //     // if (x == NONE){
    //     //     rb10_misses++;
    //     // }
    //     x = rb11.invert(y);
    //     if (x == NONE){
    //         rb11_misses++;
    //     }
    // }
    // // std::cout << "10-st table got " << rb10_misses << " misses" << std::endl;
    // std::cout << "20-nd table got " << rb11_misses << " misses" << std::endl;
    // std::cout << std::endl;
    // std::cout << std::endl;
    // std::cout << " ---- Function 5 ---- " << std::endl;
    // uint32_t p = 3371;
    // uint32_t N = 960;
    // Rainbow_Table_3 rb1(N, func_5, p, 70, 30);
    // Rainbow_Table_3 rb2(N, func_5, p, 70, 30);
    // Rainbow_Table_3 rb3(N, func_5, p, 70, 30);
    // Rainbow_Table_3 rb4(N, func_5, p, 70, 30);
    // Rainbow_Table_3 rb5(N, func_5, p, 70, 30);
    // std::unordered_set<uint32_t> image;
    // uint32_t k;
    // for (k = 0; k < N; ++k){
    //     image.insert(func_5(k));
    // }
    // std::cout << "We have " << image.size() << " elements in the image" << std::endl;
    // int rb1_misses = 0; 
    // int rb2_misses = 0; 
    // int rb3_misses = 0; 
    // int rb4_misses = 0; 
    // int rb5_misses = 0;
    // int tot_misses = 0; 
    // uint64_t x;
    // for (uint32_t y : image){
    //     bool preimage1 = true;
    //     bool preimage2 = true;
    //     bool preimage3 = true;
    //     bool preimage4 = true;
    //     bool preimage5 = true;
    //     x = rb1.invert(y);
    //     if (x == NONE){
    //         rb1_misses++;
    //         preimage1 = false;
    //     }
    //     x = rb2.invert(y);
    //     if (x == NONE){
    //         rb2_misses++;
    //         preimage2 = false;
    //     }
    //     x = rb3.invert(y);
    //     if (x == NONE){
    //         rb3_misses++;
    //         preimage3 = false;
    //     }
    //     x = rb4.invert(y);
    //     if (x == NONE){
    //         rb4_misses++;
    //         preimage4 = false;
    //     }
    //     x = rb5.invert(y);
    //     if (x == NONE){
    //         rb5_misses++;
    //         preimage5 = false;
    //     }
    //     if (!preimage1 && !preimage2 && !preimage3 && !preimage4 && !preimage5){
    //         tot_misses++;
    //     }
    // }
    // std::cout << "table 1 got " << rb1_misses << " misses" << std::endl;
    // std::cout << "table 2 got " << rb2_misses << " misses" << std::endl;
    // std::cout << "table 3 got " << rb3_misses << " misses" << std::endl;
    // std::cout << "table 3 got " << rb4_misses << " misses" << std::endl;
    // std::cout << "table 3 got " << rb5_misses << " misses" << std::endl;
    // std::cout << "Total number of misses " << tot_misses << " misses" << std::endl;
    // std::cout << std::endl;
    // std::cout << std::endl;
    



    // text length is 333614
    // window is 40001

    std::cout << " ---- Function 6 ---- " << std::endl;
    uint32_t N = func_6_prep() - 40001;
    uint64_t p = find_good_prime_for_rb_tables(N);  // 325697
    uint32_t log_N = std::log(N);
    uint32_t cbrt_N = std::cbrt(N);
    std::cout << "N = " << N << std::endl;
    std::cout << "p = " << p << std::endl;
    std::cout << "ln(N) = " << log_N << std::endl;
    uint32_t nb_chains = 2 * cbrt_N * cbrt_N;
    uint32_t chain_len = 2 * cbrt_N;
    std::cout << "nb_chains = " << nb_chains << std::endl;
    std::cout << "chain_len = " << chain_len << std::endl;
    // std::cout << "Mega RB Table completeness : " << megarb.check_completeness() << std::endl;
    Rainbow_Table_3 rb1(N, func_6, p, nb_chains, chain_len);
    Rainbow_Table_3 rb2(N, func_6, p, nb_chains, chain_len);
    Rainbow_Table_3 rb3(N, func_6, p, nb_chains, chain_len);
    Rainbow_Table_3 rb4(N, func_6, p, nb_chains, chain_len);
    Rainbow_Table_3 rb5(N, func_6, p, nb_chains, chain_len);
    Rainbow_Table_3 rb6(N, func_6, p, nb_chains, chain_len);
    Rainbow_Table_3 rb7(N, func_6, p, nb_chains, chain_len);
    Rainbow_Table_3 rb8(N, func_6, p, nb_chains, chain_len);
    Rainbow_Table_3 rb9(N, func_6, p, nb_chains, chain_len);
    Rainbow_Table_3 rb10(N, func_6, p, nb_chains, chain_len);
    std::unordered_set<uint128_t> image;
    uint32_t k;
    for (k = 0; k < N; ++k){
        image.insert(func_6(k));
    }
    std::cout << "We have " << image.size() << " elements in the image" << std::endl;
    int rb1_misses = 0; 
    int rb2_misses = 0; 
    int rb3_misses = 0; 
    int rb4_misses = 0; 
    int rb5_misses = 0;
    int rb6_misses = 0;
    int rb7_misses = 0;
    int rb8_misses = 0;
    int rb9_misses = 0;
    int rb10_misses = 0;
    int tot_misses = 0; 
    uint64_t x;
    int u = 0;
    for (uint128_t y : image){
        u++;
        // std::cout << "u = " << u << std::endl;
        if (u%10000 == 0){
            std::cout << "u = " << u << std::endl;
            // std::cout << "table 1 got " << rb1_misses << " misses" << std::endl;
            // std::cout << "table 2 got " << rb2_misses << " misses" << std::endl;
            // std::cout << "table 3 got " << rb3_misses << " misses" << std::endl;
            // std::cout << "table 4 got " << rb4_misses << " misses" << std::endl;
            // std::cout << "table 5 got " << rb5_misses << " misses" << std::endl;
            // std::cout << "table 6 got " << rb6_misses << " misses" << std::endl;
            // std::cout << "table 7 got " << rb7_misses << " misses" << std::endl;
            // std::cout << "table 8 got " << rb8_misses << " misses" << std::endl;
            // std::cout << "table 9 got " << rb9_misses << " misses" << std::endl;
            // std::cout << "table 10 got " << rb10_misses << " misses" << std::endl;
            std::cout << "Total number of misses " << tot_misses << " misses" << std::endl;
        } 
        bool preimage1 = true;
        bool preimage2 = true;
        bool preimage3 = true;
        bool preimage4 = true;
        bool preimage5 = true;
        bool preimage6 = true;
        bool preimage7 = true;
        bool preimage8 = true;
        bool preimage9 = true;
        bool preimage10 = true;
        x = rb1.invert(y);
        if (x == NONE){
            rb1_misses++;
            preimage1 = false;
        }
        x = rb2.invert(y);
        if (x == NONE){
            rb2_misses++;
            preimage2 = false;
        }
        x = rb3.invert(y);
        if (x == NONE){
            rb3_misses++;
            preimage3 = false;
        }
        x = rb4.invert(y);
        if (x == NONE){
            rb4_misses++;
            preimage4 = false;
        }
        x = rb5.invert(y);
        if (x == NONE){
            rb5_misses++;
            preimage5 = false;
        }
        x = rb6.invert(y);
        if (x == NONE){
            rb6_misses++;
            preimage6 = false;
        }
        x = rb7.invert(y);
        if (x == NONE){
            rb7_misses++;
            preimage7 = false;
        }
        x = rb8.invert(y);
        if (x == NONE){
            rb8_misses++;
            preimage8 = false;
        }
        x = rb9.invert(y);
        if (x == NONE){
            rb9_misses++;
            preimage9 = false;
        }
        x = rb10.invert(y);
        if (x == NONE){
            rb10_misses++;
            preimage10 = false;
        }
        if (!preimage1 && !preimage2 && !preimage3 && !preimage4 && !preimage5 && !preimage6&& !preimage7&& !preimage8&& !preimage9&& !preimage10){
            tot_misses++;
        }
    }
    std::cout << "u = " << u << std::endl;
    std::cout << "table 1 got " << rb1_misses << " misses" << std::endl;
    std::cout << "table 2 got " << rb2_misses << " misses" << std::endl;
    std::cout << "table 3 got " << rb3_misses << " misses" << std::endl;
    std::cout << "table 4 got " << rb4_misses << " misses" << std::endl;
    std::cout << "table 5 got " << rb5_misses << " misses" << std::endl;
    std::cout << "table 6 got " << rb6_misses << " misses" << std::endl;
    std::cout << "table 7 got " << rb7_misses << " misses" << std::endl;
    std::cout << "table 8 got " << rb8_misses << " misses" << std::endl;
    std::cout << "table 9 got " << rb9_misses << " misses" << std::endl;
    std::cout << "table 10 got " << rb10_misses << " misses" << std::endl;
    std::cout << "Total number of misses " << tot_misses << " misses" << std::endl;
    std::cout << std::endl;
    std::cout << std::endl;
    return 0;
}

// Old run for function 6 with text of size ~ 160_000:
// Output of function 6: 
// table 1 got 15411 misses
// table 2 got 14940 misses
// table 3 got 15231 misses
// table 4 got 15024 misses
// table 5 got 15310 misses
// table 6 got 15049 misses
// Total number of misses 0 misses

// For some purposes, when the number of misses is low, even 3 tables are sufficient.
// Through multiple tests with 6 rb tables, it has occurred that 1 miss is still possible.
// Even though the image is typically ~133_000 (depending on the base chosen randomly)
// We strongly avoid false-negatives.