#pragma once
#include <iostream>
#include <string>
#include <vector>
#include <stack>
#include <array>
#include <memory>


struct rads_simple_1{
    size_t tau;
    size_t buffer_start;
    size_t buffer_size;
    size_t l_match;
    size_t i_match;
    const std::string_view R;
    const std::string_view T;
    rads_simple_1(std::string& ref, std::string& text, size_t tau) : R {ref}, T {text}, tau {tau},
    l_match {0}, i_match {0}, buffer_start {0}, buffer_size{0} {
        if (T.size() == 0){
            throw std::runtime_error("Text is empty --> undefined behaviour. Initialization ceased");
        }
    }
    char text_at(size_t j){
        if (j < l_match){
            return R[i_match+j];
        } else if (j-l_match < buffer_size){
            return T[buffer_start + j-l_match];
        } else{
            throw std::runtime_error("Bad indexing ... Access denied");   // Remove if the code is correct
        }
    }
    char rev_text_at(size_t j){
        if (j < buffer_size){
            return T[buffer_start+buffer_size-1-j];
        } else if (j-buffer_size < l_match){
            return R[i_match+l_match-1-(j-buffer_size)];
        } else{
            throw std::runtime_error("Bad indexing ... Access denied");  // Remove if the code is correct
        }
    }
    char ref_at(size_t j){
        if (j < R.size()){
            return R[j];
        } else{
            throw std::runtime_error("Bad indexing ... Access denied");  // Remove if the code is correct
        } 
    }
    char rev_ref_at(size_t j){
        if (j < R.size()){
            return R[R.size()-1-j];
        } else{
            throw std::runtime_error("Bad indexing ... Access denied");   // Remove if the code is correct
        }
    }
    std::tuple<size_t, size_t, size_t> UpdateMS_rev(size_t l, size_t s, size_t p){
        if (l == 0){
            return {1,0,1};
        }
        size_t i = l;
        while (i < l + 1){
            if (rev_text_at(i-p) < rev_text_at(i)){               
                i = i - (i>=s ? (i-s)%p : p-((s-i)%p));  
                s = i;                                   
                p = 1;                                   
            } else if(rev_text_at(i-p) > rev_text_at(i)){
                p = i - s + 1;
            }
            i++;
        }
        return {l+1, s, p};
    }
    std::pair<size_t, size_t> longest_suffix_match(size_t query_st = 0, size_t query_len = 0){
        size_t i = 0, p = 0, s = 0, l = 0, l_max = 0;
        size_t m = query_len == 0 ? l_match + buffer_size : query_len;
        size_t n = R.size(); 
        size_t ans_pos = 0;
        if (m == 0){
            return std::make_pair(0, 0);
        }
        while (i < n){
            while (i+l<n && l<m && rev_ref_at(i+l)==rev_text_at(l)){
                std::tie(l, s, p) = UpdateMS_rev(l, s, p);
            }
            if (l >= l_max){
                ans_pos = i;
                l_max = l;
            }
            bool s_match = true;
            for (size_t b = 0; b < s; ++b){ 
                if(rev_text_at(b)!=rev_text_at(p+b)){
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
    void update(){
        if (buffer_size < tau && buffer_start+buffer_size < T.size()){
            buffer_size++;
            return ;
        }
        std::tie(i_match, l_match) = longest_suffix_match();
        i_match = R.size() - i_match - l_match;
        buffer_start += buffer_size;
        buffer_size = std::min(T.size()-buffer_start, (size_t) 1);
    }
    void fill(){  // run only once at the start, all subsequent fillings are made by invoking update()
        buffer_size = std::min(T.size()-buffer_start,tau);
    }
    void print(){
        std::cout << "(i_match, l_match) = " << "(" << i_match << " , " 
        << l_match << ")" << std::endl;
        std::cout << "--Match in R--" << std::endl;
        std::cout << R.substr(i_match, l_match) << std::endl;
        std::cout << "buffer_start = " << buffer_start << std::endl;
        std::cout << "buffer_size = " << buffer_size << std::endl;
        std::cout << "--T supported--" << std::endl;
        std::cout << T.substr(buffer_start, buffer_size) << std::endl;
    }
};


struct rads_simple_2{
    const std::string_view R;
    const std::string_view T;
    size_t tau;
    size_t l_match;
    size_t i_match;
    size_t buffer_start;
    size_t buffer_size;
    std::vector<char> current_buff;
    size_t N; // Currently N = T.size(), but it'll become the global size of the text data, whereas T becomes the buffer for RIO-reading of the data.
    rads_simple_2(std::string& ref, std::string& text, size_t tau) : R {ref}, T {text}, tau {tau},
    l_match {0}, i_match {0}, buffer_start {0}, buffer_size{0}, current_buff(tau), N{T.size()} {
        if (T.size() == 0){
            throw std::runtime_error("Text is empty --> undefined behaviour. Initialization ceased");
        }

    }
    //void update_text();
    char text_at(size_t j){
        if (j < l_match){
            return R[i_match+j];
        } else if (j-l_match < buffer_size){
            return current_buff[j-l_match];
        } else{
            throw std::runtime_error("Bad indexing ... Access denied");   // Remove if the code is correct
        }
    }
    char rev_text_at(size_t j){
        if (j < buffer_size){
            return current_buff[buffer_size-1-j];
        } else if (j-buffer_size < l_match){
            return R[i_match+l_match-1-(j-buffer_size)];
        } else{
            throw std::runtime_error("Bad indexing ... Access denied");  // Remove if the code is correct
        }
    }
    char ref_at(size_t j){
        if (j < R.size()){
            return R[j];
        } else{
            throw std::runtime_error("Bad indexing ... Access denied");  // Remove if the code is correct
        } 
    }
    char rev_ref_at(size_t j){
        if (j < R.size()){
            return R[R.size()-1-j];
        } else{
            throw std::runtime_error("Bad indexing ... Access denied");   // Remove if the code is correct
        }
    }
    std::tuple<size_t, size_t, size_t> UpdateMS_rev(size_t l, size_t s, size_t p){
        if (l == 0){
            return {1,0,1};
        }
        size_t i = l;
        while (i < l + 1){
            if (rev_text_at(i-p) < rev_text_at(i)){               
                i = i - (i>=s ? (i-s)%p : p-((s-i)%p));  
                s = i;                                   
                p = 1;                                   
            } else if(rev_text_at(i-p) > rev_text_at(i)){
                p = i - s + 1;
            }
            i++;
        }
        return {l+1, s, p};
    }
    std::pair<size_t, size_t> longest_suffix_match(size_t query_st = 0, size_t query_len = 0){
        size_t i = 0, p = 0, s = 0, l = 0, l_max = 0;
        size_t m = query_len == 0 ? l_match + buffer_size : query_len;
        size_t n = R.size(); 
        size_t ans_pos = 0;
        if (m == 0){
            return std::make_pair(0, 0);
        }
        while (i < n){
            while (i+l<n && l<m && rev_ref_at(i+l)==rev_text_at(l)){
                std::tie(l, s, p) = UpdateMS_rev(l, s, p);
            }
            if (l >= l_max){
                ans_pos = i;
                l_max = l;
            }
            bool s_match = true;
            for (size_t b = 0; b < s; ++b){ 
                if(rev_text_at(b)!=rev_text_at(p+b)){
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
    char receive(){
        if (buffer_start+buffer_size < T.size()){
            char c = T[buffer_start+buffer_size];
            return c;
        }
        // } else{
        //     if (buffer_start+buffer_size < N){
        //         update_text();
        //         char c = T[0];
        //         buffer_start = 0;
        //         return c;
        //     }
        //     return;
        // }
    }
    void update(){
        std::tie(i_match, l_match) = longest_suffix_match();
        i_match = R.size() - i_match - l_match;
        buffer_start += buffer_size;
        buffer_size = 0;
        while (buffer_size < tau && buffer_start+buffer_size < N){
            current_buff[buffer_size] = receive();
            buffer_size++;
        }
    }
    void fill(){
        size_t target = std::min(T.size() - buffer_start, tau);
        while (buffer_size < target){
            current_buff[buffer_size] = receive();
            buffer_size++;
        }
    }
    void print(){
        std::cout << "(i_match, l_match) = " << "(" << i_match << " , " 
        << l_match << ")" << std::endl;
        std::cout << "--Match in R--" << std::endl;
        std::cout << R.substr(i_match, l_match) << std::endl;
        std::cout << "buffer_start = " << buffer_start << std::endl;
        std::cout << "buffer_size = " << buffer_size << std::endl;
        std::cout << "--T supported--" << std::endl;
        for (auto s : current_buff){
            std::cout << s;
        }
        std::cout << std::endl;
    }
};
 