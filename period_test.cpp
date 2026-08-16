#include <iostream>
#include "compute_period.hpp"

int main(){
    std::string w1 = "abbabbabbabbab";
    access_func access1 = [&w1](uint32_t k)->unsigned char {
        return (unsigned char)w1[k];
    };
    std::string w2 = "baababaabaababaabaababaabaababaa";
    access_func access2 = [&w2](uint32_t k)->unsigned char {
        return (unsigned char)w2[k];
    };
    std::string w3 = "abcdefghkelmnufh";
    access_func access3 = [&w3](uint32_t k)->unsigned char {
        return (unsigned char)w3[k];
    };
    std::cout << "Period of w1 is " << compute_period(0, w1.size(), access1) << "\nReal period is equal to 3" << std::endl;
    std::cout << "Period of w2 is " << compute_period(0, w2.size(), access2) << "\nReal period is equal to 8" << std::endl;
    std::cout << "Period of w3 is " << compute_period(0, w3.size(), access3) << "\nReal period is equal to 16" << std::endl;
    return 0;
}