#include <iostream>
#include <sdsl/suffix_trees.hpp>
#include <vector>
#include <stack>
#include <utility>
#include "pm&rlz_utils.hpp"

// TODO :
// Add std::move logic, because currently std::strings 
// are copied each time PM is invoked

#define SIGMA 2

constexpr char ALPHABET[] = {'a', 'b'};

using table_t = std::vector<std::vector<size_t>>;

using node_type = sdsl::cst_sada<>::node_type;

std::vector<size_t> calculate_table_of_deepest_nodes(sdsl::cst_sada<>& cst, table_t& marked)
{
    std::cout << "deepest nodes reached" << std::endl;
    std::vector<size_t> table_of_best(cst.csa.size());
    std::stack<std::tuple<node_type, size_t, node_type>> stack;
    node_type root = cst.root();
    stack.push({root, 0, root});
    while(!stack.empty()){
        std::tuple state = stack.top();
        node_type v = std::get<0>(state);
        size_t val = std::get<1>(state);
        node_type candidate = std::get<2>(state);
        stack.pop();
        size_t id_v = cst.id(v);
        if (marked[id_v][0] == 1){
            size_t tmp = marked[id_v][2] + marked[id_v][3];
            if (tmp > val){
                candidate = v;
                val = tmp;
            }
        }
        if (cst.is_leaf(v)){
            size_t suff_numb = cst.sn(v);
            table_of_best[suff_numb] = cst.id(candidate);
        } else {
            for (uint8_t y = 0; y < SIGMA; y++){
                node_type w = cst.child(v, ALPHABET[y]);
                if (w != root){   
                    stack.push({w, val, candidate});
                }
            }
        }
    }

    return table_of_best;
}


void write_into_marked(size_t& id_v, size_t& v_depth, size_t& i, size_t cur_match, table_t& marked){
    if (cur_match >= marked[id_v][2]){
        marked[id_v][0] = 1;
        marked[id_v][1] = i; // The starting position of the suffix of R concerned rn
        marked[id_v][2] = cur_match; // Offset on the edge = number of matched chars on it
        marked[id_v][3] = v_depth;
    }
}


void extend(node_type& v, size_t& i, table_t& marked, sdsl::cst_sada<>& suff_tree, size_t& m){
    while (i < m){
        size_t id_v = suff_tree.id(v);
        size_t v_depth = suff_tree.depth(v);
        write_into_marked(id_v, v_depth, i, 0, marked);
        i++;
        v = suff_tree.sl(v);
        std::cout << "i = " << i << std::endl;
    }
    return ;
};


std::vector<std::tuple<size_t, size_t, size_t>> rlz_compress(const std::string& text, const std::string& ref, size_t tau)
{
    size_t m = ref.size();
    size_t n = text.size();
    size_t global_offset = 0;
    std::vector<std::tuple<size_t, size_t, size_t>> factors;
    size_t b = 0;
    size_t block_size;
    while (b < n){
        std::cout << "b = " << b << std::endl;
        block_size = std::min(tau,n-b);
        std::string block = text.substr(b, block_size);
        sdsl::cst_sada<> suff_tree;
        sdsl::construct_im(suff_tree, block, 1);
        std::cout << block << std::endl;
        std::cout<< "Suff tree construction is done" <<std::endl;
        table_t marked(suff_tree.nodes(), std::vector<size_t> (4));
        size_t i = 0;
        size_t j = 0;
        node_type v = suff_tree.root();
        while (i < m){
            while (j < m){
                std::cout<< "i = "<< i << "     j = " << j <<std::endl;
                if (suff_tree.is_leaf(v)){
                    std::cout<< "Reached a leaf" << std::endl;
                    size_t id_v = suff_tree.id(v);
                    size_t v_depth = suff_tree.depth(v);
                    write_into_marked(id_v, v_depth, i, 0, marked);
                    v = suff_tree.sl(v);
                    i++;
                    j = i;
                    break;
                }
                // Trying to read a new char
                char c = ref[j];
                std::cout<< "Trying to read a new char = " << c << std::endl;
                node_type w = suff_tree.child(v, c); 
                size_t cur_match = 0;
                size_t v_depth = suff_tree.depth(v);
                int edge_length = suff_tree.depth(w) - v_depth;
                if (edge_length <= 0){
                    size_t id_v = suff_tree.id(v);
                    std::cout<< "No match on the edge i = "<< i << "     j = " << j <<std::endl;
                    write_into_marked(id_v, v_depth, i, 0, marked);
                    j = std::max(i + marked[suff_tree.id(v)][3], i + 1);
                    ++i;
                    v = suff_tree.sl(v);
                    if (j >= m) 
                        extend(v, i, marked, suff_tree, m);
                    break;
                }
                while (cur_match < edge_length && j + cur_match < m && suff_tree.edge(w,v_depth + cur_match + 1) == ref[j + cur_match]){  //1-based indexing on edges
                    ++cur_match;
                    std::cout<< "cur_match = " << cur_match << std::endl;
                }
                if (cur_match < edge_length){
                    std::cout<< "Partial match on the edge i = "<< i << "     j = " << j <<std::endl;
                    size_t id_v = suff_tree.id(v);
                    write_into_marked(id_v, v_depth, i, cur_match, marked);
                    j = std::max(i + marked[suff_tree.id(v)][3], i + 1);
                    ++i;
                    v = suff_tree.sl(v);
                    break;
                } else{
                    std::cout<< "Complete match on the edge" << std::endl;
                    j+= edge_length;
                    std::cout<< " i = " << i << "   j = "<< j << std::endl;
                    v = w;
                    if (j >= m){
                        extend(v, i, marked, suff_tree, m);
                        std::cout<< " i = " << i << "   j = "<< j << std::endl;
                        break;
                    }
                }
            }
        }
        std::cout<<"Going to PM" << std::endl;
        std::vector<size_t> table_of_best = calculate_table_of_deepest_nodes(suff_tree, marked);
        std::cout <<"marked[table_of_best[0]][2] + marked[table_of_best[0]][3] =" << marked[table_of_best[0]][2] + marked[table_of_best[0]][3] << std::endl;
        int p = 0; 
        while (p + global_offset < block.size()){
            // std::cout<< "b + p + global_offset = " << b+p+global_offset << std::endl;
            size_t idv = table_of_best[p+global_offset];
            size_t length_of_factor = marked[idv][2] + marked[idv][3];
            size_t start_of_factor = marked[idv][1];
            factors.push_back(std::make_tuple(p+global_offset,start_of_factor, length_of_factor));
            p += length_of_factor;
            // std::cout<< "First factor : (pos, len) = " << std::get<1>(factors.front()) << " " << std::get<2>(factors.front()) << std::endl;
            // std::cout<< "New factor : (pos, len) = " << start_of_factor << " " << length_of_factor << std::endl;
        }
        // Here we assume that have found at least one factor. (Meaning that forall c\in ALPHABET, c is contained in R)
        size_t x = b + std::get<0>(factors.back());
        while (b < n && x + std::get<2>(factors.back()) == (b+block_size)){
            std::cout<< "PM" << std::endl;
            b += block_size;
            block_size = std::min(n-b, tau);
            std::cout << "Within PM: block_size = " << block_size << std::endl;
            block = text.substr(b, block_size);
            pm_rlz PM{std::get<1>(factors.back()), std::get<2>(factors.back()), ref, block};
            std::cout << "PM is constructed" << std::endl;
            auto pair = PM.longest_prefix_match();
            std::cout << "PM is done" << std::endl;
            std::get<1>(factors.back()) = pair.first;
            std::get<2>(factors.back()) = pair.second;
        }
        global_offset = x + std::get<2>(factors.back()) - b;
    }
    std::cout << std::endl;
    std::cout << std::endl;
    std::cout << std::endl;
    std::cout << "--- FACTORS ---"<< std::endl;
    for (auto s: factors){
        std::cout<< "Start in R =" << std::get<1>(s) <<
        " length = " << std::get<2>(s) << std::endl;
    }
    return factors;
}

// Resulting complexity O(m + tau*log(tau)) = O(m), if tau is chosen properly

int main(int argc, char* argv[]) {
    size_t tau = 4;
    std::string ref = "abba";
    std::string text = "abba";
    std::vector<std::tuple<size_t, size_t, size_t>> factorization = rlz_compress(text, ref, tau);
    return 0;
}
