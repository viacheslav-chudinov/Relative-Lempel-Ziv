#include "rads_ds.hpp"

void rads_in_action(std::string& T, std::string& R,size_t tau){
    size_t n = T.size();
    rads_primitive rads(R, T, tau);
    while (rads.buffer_size > 0){
        std::cout << std::endl;
        std::cout << " -- Current state is -- " << std::endl;
        std::cout << "Match in R:" << std::endl;
        std::cout << R.substr(rads.i_match, rads.l_match) << std::endl;;
        std::cout << "Current buffer's state is" << std::endl;
        for (auto s = rads.buffer_start; s < rads.buffer_start + rads.buffer_size; ++s){
            std::cout << T[s];
        } 
        std::cout << std::endl;
        rads.update();
    }
    std::cout << "At the end : \n Current buffer's state is" << std::endl;
    for (auto s = rads.buffer_start; s < rads.buffer_start + rads.buffer_size; ++s){
        std::cout << T[s];
    }
    std::cout << std::endl;
    std::cout << "rads.buffer_size = " << rads.buffer_size << std::endl;
}

int main(){
    std::string T = "babbbbaabbab";
    std::string R = "bbaabab";
    size_t tau = 3;
    rads_in_action(T, R, tau);
    return 0;
}