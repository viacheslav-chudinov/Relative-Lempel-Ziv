//rlz_test+primitive_rlz.cpp
#include <iostream>
#include "pm&rlz_utils.hpp"
#include <random>


factors_t rlz_compress_primitive(std::string& text, std::string& ref)
{
    size_t m = ref.size();
    size_t n = text.size();
    factors_t factors;
    size_t b = 0;
    while (b < n) {
        size_t best_i = 0;
        size_t best_len = 0;
        for (size_t i = 0; i < m; ++i) {
            size_t len = 0;
            while (i + len < m && b + len < n && text[b + len] == ref[i + len]) {
                ++len;
            }
            if (len > best_len) {
                best_len = len;
                best_i = i;
            }
        }
        factors.push_back(std::make_tuple(-1, best_i, best_len));
        b += best_len;
    }
    return factors;
}



int main(){
    std::vector<unsigned > nb_tests{2, 3, 3};
    std::vector<unsigned > ref_lengths{16, 25, 36};
    std::vector<unsigned> taus{4, 5, 6};
    unsigned int text_length = 50;
    for (int Q = 0; Q < 2; ++Q){
        unsigned int nb_of_tests = nb_tests[Q];
        unsigned int ref_length = ref_lengths[Q];
        size_t tau = taus[Q];
        for (std::size_t n = 0; n < nb_of_tests; ++n)
        {   
            std::cout << "Ref length is  " << ref_length << "\n And test number " << n << "started" << std::endl;
            std::mt19937 generator(std::random_device{}());
            std::uniform_int_distribution<std::size_t> distribution(0, SIGMA - 2);
            std::string ref = "";
            std::string text = "";
            for (unsigned int i = 0; i < ref_length; ++i){
                std::size_t number = distribution(generator);
                ref += ALPHABET[number];
            }
            for (unsigned int i = 0; i < text_length; ++i){
                std::size_t number = distribution(generator);
                text += ALPHABET[number];
            }
            erlz_simple RLZ_simple{tau, ref, text};
            factors_t factors = RLZ_simple.rlz_compress();
            factors_t true_factors = rlz_compress_primitive(text, ref);
            size_t l1, l2;
            l1 = factors.size();
            l2 = true_factors.size();
            if (l1 != l2){
                std::cout << "WRONG FACTORIZATION CAUGHT" << std::endl;
                throw("WRONG FACTORIZATION CAUGHT");
            }            
            for (size_t u = 0; u < l1; ++u){
                auto factor = factors[u];
                auto true_factor = true_factors[u];
                size_t lenf1, lenf2;
                lenf1 = std::get<2>(factor);
                lenf2 = std::get<2>(true_factor);
                size_t st1, st2;
                st1 = std::get<1>(factor);
                st2 = std::get<1>(true_factor);
                if (lenf1 != lenf2){
                    std::cout << "WRONG FACTORIZATION CAUGHT" << std::endl;
                    throw("WRONG FACTORIZATION CAUGHT");
                }
                for (size_t v = 0; v < lenf1; ++v){
                    if (ref[st1 + v] != ref[st2+v]){
                        std::cout << "WRONG FACTORIZATION CAUGHT" << std::endl;
                        throw("WRONG FACTORIZATION CAUGHT");
                    }
                }
            }
            std::cout << "Ref length is  " << ref_length << "\n And test number n has been completed" << std::endl;

        }
    }
    std::cout << "OK  -- RLZ_SIMPLE IS CORRECT" << std::endl;

    std::vector<std::string> refs{"abbedced", "bbbedabaddeedbdccdace", "aabbabaabbbababbaaaaaaabbabaabbabbbb"};
    std::vector<std::string> texts = {"bbebedcedabbeedcbcebedcedab", "edabaddeedbdccdedbdbbeeeacdabadccdacbbbaddeedbdccdace", "baabbbabababbaaaaaaabbaabbabaabbabbabbabaabbbababbababbaaaaaaabbbabaabbabbbbbaaaaaaabbabaabaabbbababbaaaaaaabbabaa"};
    std::string ref = "abbedced";
    std::string text = "bbedceabbedcedadc";
    std::vector<size_t> tausn = {3, 5, 10};
    for (int i =0 ; i < 2; ++i){
        std::string ref = refs[i];
        std::string text = texts[i];
        size_t tau = tausn[i];
        erlz_simple RLZ_simple{tau, ref, text};
        factors_t factors = RLZ_simple.rlz_compress();
        factors_t true_factors = rlz_compress_primitive(text, ref);
        size_t l1, l2;
        l1 = factors.size();
        l2 = true_factors.size();
        if (l1 != l2){
            std::cout << "WRONG FACTORIZATION CAUGHT" << std::endl;
            throw("WRONG FACTORIZATION CAUGHT");
        }            
        for (size_t u = 0; u < l1; ++u){
            auto factor = factors[u];
            auto true_factor = true_factors[u];
            size_t lenf1, lenf2;
            lenf1 = std::get<2>(factor);
            lenf2 = std::get<2>(true_factor);
            size_t st1, st2;
            st1 = std::get<1>(factor);
            st2 = std::get<1>(true_factor);
            if (lenf1 != lenf2){
                std::cout << "WRONG FACTORIZATION CAUGHT" << std::endl;
                throw("WRONG FACTORIZATION CAUGHT");
            }
            for (size_t v = 0; v < lenf1; ++v){
                if (ref[st1 + v] != ref[st2+v]){
                    std::cout << "WRONG FACTORIZATION CAUGHT" << std::endl;
                    throw("WRONG FACTORIZATION CAUGHT");
                }
            }
        }
        std::cout << "OK  -- RLZ_SIMPLE IS CORRECT" << std::endl;
    }
    return 0;
}

