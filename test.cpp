#include <iostream>
#include <vector>
#include <cstddef>
#include <utility>
#include <string> 

// 1. Создать папку сборки
// mkdir build
// cd build
// 2. Сгенерировать проект
// cmake ..
// 3. Собрать
// make

using state_t = std::tuple<size_t, size_t, size_t>;

state_t UpdateMS(std::string& Y, size_t l, size_t s, size_t p){
    if (l == 0){
        return {1,0,1};
    }
    size_t i = l;
    while (i < l + 1){
        if (Y[i-p] < Y[i]){
            i = i - (i>=s ? (i-s)%p : p-((s-i)%p));
            s = i;
            p = 1;
        } else if(Y[i-p] > Y[i]){
            p = i - s + 1;
        }
        i++;
    }
    return {l+1, s, p};
}
std::pair<size_t, size_t> longest_prefix_match(std::string& R, std::string T){
    size_t i = 0, p = 0, s = 0, l = 0, l_max = 0;
    size_t m = R.size();
    size_t n = T.size(); 
    size_t ans_pos = 0;
    while (i < n){
        while (i+l<n && l<m && T[i+l]==R[l]){
            std::tie(l, s, p) = UpdateMS(R, l, s, p);
        }
        if (l >= l_max){
            ans_pos = i;
            l_max = l;
        }
        bool s_match = true;
        for (size_t b = 0; b < s; ++b){ 
            if(R[b]!=R[p+b]){
                s_match = false;
                break;                    
            }
        }
        if (p > 0 && p <= l/3 && s_match){
            i += p;
            l -= p;
        } else{
            i += l/3 + 1;
            l = 0; s = 0; p = 0;
        }
        }
        return std::make_pair(ans_pos, l_max);
    }
int main(){
    std::string R = "abbcba";
    std::string T = "cabdabeabbcgghdabbcbbbabbay";
    std::pair res = longest_prefix_match(R, T);
    std::cout << " position in T found = " << res.first << std::endl;
    std::cout << " length of the longest occurrence = " << res.second << std::endl;
    return 0;
}