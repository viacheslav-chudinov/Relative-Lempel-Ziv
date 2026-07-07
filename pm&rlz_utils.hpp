#pragma once
#include <iostream>
#include <sdsl/suffix_trees.hpp>
#include <vector>
#include <utility>


using factors_t = std::vector<std::tuple<size_t, size_t, size_t>>;

inline constexpr size_t SIGMA = 6;
inline char ALPHABET[SIGMA] = {'a','b','c','d','e','\0'};



// Algorithm taken from https://www.stringology.org/event/2013/p15.html
struct pm_rlz{
    size_t start;
    size_t length;
    std::string& R;
    std::string& block;
    char ref_at(size_t pos){
        if (pos >= length && pos < length + block.size()){
            return block[pos - length];
        } else if (pos < length && start+pos<R.size()){
            return R[start+pos];
        } else{
            std::cout << "start+pos = " << start+pos << std::endl;
            std::cout << "length = " << length << std::endl;
            throw std::runtime_error("ref_at bad indexing");
        }
    }
    char txt_at(size_t pos){
        if (pos<R.size()){
            return R[pos];
        } else{
            std::cout << "pos = " << pos <<std::endl; 
            throw std::runtime_error("text_at bad indexing");
        }
    }
    std::tuple<size_t, size_t, size_t> UpdateMS(size_t l, size_t s, size_t p){
        if (l == 0){
            return {1,0,1};
        }
        size_t i = l;
        while (i < l + 1){
            if (ref_at(i-p) < ref_at(i)){
                i = i - (i>=s ? (i-s)%p : p-((s-i)%p));
                s = i;
                p = 1;
            } else if(ref_at(i-p) > ref_at(i)){
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
    
}; 

using table_t = std::vector<std::vector<size_t>>;

using node_type = sdsl::cst_sada<>::node_type;

struct erlz_simple{
    size_t tau;
    std::string& ref;
    std::string& text;

    std::vector<size_t> calculate_table_of_deepest_nodes(sdsl::cst_sada<>& cst, table_t& marked);

    void write_into_marked_1(size_t& id_v, size_t& v_depth, size_t& i, table_t& marked);
    void write_into_marked_0(size_t& id_w, size_t& v_depth, size_t& i, size_t cur_match, table_t& marked);
    void extend(node_type& v, size_t& i, table_t& marked, sdsl::cst_sada<>& suff_tree, size_t& m);

    factors_t rlz_compress();
};