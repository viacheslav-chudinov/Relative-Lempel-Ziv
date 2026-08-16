#include <iostream>
#include "primes.hpp"

int main(){
    std::vector<uint64_t> primes = {37, 2789, 3041, 17581, 18773, 21787, 
    21817, 118463, 136069, 391711, 431857, 497867, 529313, 608389};  // prime numbers
    std::vector<uint64_t> composites = {35, 2793, 2805, 3479, 18453, 21931,
    29519, 119729, 348037, 431261, 520417, 529093};
    for (auto prime : primes){
        bool flag = is_prime(prime);
        if (!flag){
            throw std::runtime_error("Prime " + std::to_string(prime) + " was identified as a non-prime");
        }
        std::cout << prime << " was correctly identified" << std::endl;
    }
    for (auto composite : composites){
        bool flag = is_prime(composite);
        if (flag){
            throw std::runtime_error("Composite " + std::to_string(composite) + " was declared to be a prime");
        }
        std::cout << composite << " was correctly identified" << std::endl;
    }

    return 0;
}