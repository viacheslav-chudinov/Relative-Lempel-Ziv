#pragma once
#include <cstddef>
#include <utility>
#include <string> 

using state_t = std::tuple<size_t, size_t, size_t>;

// Algorithm taken from https://www.stringology.org/event/2013/p15.html

struct pm_rlz{
    size_t start;
    size_t length;
    const std::string& R;
    const std::string& block;
    char ref_at(size_t pos){
        if (pos>=length){
            return block[pos - length];
        } else{
            return R[start+pos];
        }
    }
    char txt_at(size_t pos){
        return R[pos];
    }
    state_t UpdateMS(size_t l, size_t s, size_t p){
        if (l == 0){
            return {1,0,1};
        }
        size_t i = l;
        while (i < l + 1){
            if (ref_at(i-p) > ref_at(i)){
                i = i - ( i-s >=0 ? (i-s)%p : (s-i)%p);
                s = i;
                p = 1;
            } else if(ref_at(i-p) < ref_at(i)){
                p = i - s + 1;
            }
            i++;
        }
        return {l+1, s, p};
    }
    std::pair<size_t, size_t> longest_prefix_match(){
        size_t i = 0, p = 0, s = 0, l = 0, l_max = 0;
        size_t m = length + block.size();
        size_t n = R.size(); 
        size_t ans_pos = 0;
        while (i < n){
            while (i+l<n && l<m && txt_at(i+l)==ref_at(l)){
                std::tie(l, s, p) = UpdateMS(l, s, p);
            }
            if (l >= l_max){
                ans_pos = i;
                l_max = l;
            }
            bool s_match = true;
            for (size_t b = 0; b < s; ++b){  // Sketchy
                if(ref_at(b)!=ref_at(p+b)){
                    s_match = false;
                    break;                    
                }
            }
            if (p <= l/3 && s_match){
                i += p;
                l -= p;
            } else{
                i += l/3 + 1;
                l = 0; s = 0; p = 0;
            }
        }
        return std::make_pair(l_max, ans_pos);
    }
    
}; 

