#include <iostream>
#include <sdsl/suffix_trees.hpp>
#include <vector>
#include <utility>

#define SIGMA 26;

using table_t = std::vector<std::vector<size_t>>;

using node_type = sdsl::cst_sada<>::node_type;

std::vector<size_t> calculate_table_of_deepest_nodes(sdsl::cst_sada<>& cst, table_t& marked)
{
    std::vector<size_t> table_of_best(cst.csa.size());
    std::stack<node_type> stack;
    std::stack<std::pair<<node_type, size_t>> path_best;
    node_type root = cst.root();
    stack.push_back(root);
    while(!stack.empty()){
        node_type v = stack.top();
        stack.pop();
        size_t id_v = cst.id(v);
        size_t val = path_best.top().second;
        node_type candidate = path_best.top().first;
        path_best.pop();
        if (marked[id_v][0] == 1){
            size_t tmp = marked[id_v][2] + marked[id_v][3];
            if (tmp > val){
                candidate = v;
                val = tmp;
            }
        }
        if (cst.is_leaf(v)){
            size_t sa_index = cst.sn(v);
            table_of_best[sa_index] = cst.id(candidate);
        } else {
            for (size_t c = 0; c < SIGMA; ++c){
                node_type w = cst.child(v, c);
                if (w != root)   
                    stack.push(w);
                    path_best.push({candidate, val});
            }
        }
    }

    return table_of_best;
}

rlz_compress(const std::string& text, const std::string& ref, size_t tau)
{
    size_t m = ref.size();
    size_t n = text.size();
    size_t offset = 0;
    std::vector<std::tuple<size_t, size_t, size_t>> factors;
    size_t b = 0;
    size_t block_size;
    while (b < n){
        std::string block = text.substr(b, std::min(tau,n-b));
        block_size = block.size();
        block.push_back("$");
        sdsl::cst_sada<> suff_tree;
        sdsl::construct_im(suff_tree, block, 1);

        // TODO : FIND OUT THE WAY TO STORE THE MARKED POSITIONS SUCH THAT THEY ARE EASILY ACCESSIBLE DURING LAST DFS TRAVERSAL
        table_t marked(suff_tree.nodes(), std::vector<size_t> (4));
        size_t i = 0;
        size_t j = i;
        node_type v = suff_tree.root();

        while (i < m){
            while (j < m){
                if (suffix_trees.is_leaf(v)){
                    marked[suff_tree.id(v)][0] = 1; // 1 = being marked as deepest on a go
                    marked[suff_tree.id(v)][1] = i; // The starting position of the suffix of R concerned rn
                    marked[suff_tree.id(v)][2] = 0; // Offset on the edge = number of matched chars on it
                    marked[suff_tree.id(v)][3] = suff_tree.depth(v);
                    v = suff_tree.sl(v);
                    i++;
                    j = i;
                    break;
                }
                // Trying to read a new char
                char c = ref[j];
                node_type w = suff_tree.child(v, c); 
                int cur_match = 0;
                size_t v_depth = suff_tree.depth(v);
                int edge_length = suff_tree.depth(w) - v_depth;
                while (cur_match < edge_length && suff_tree.edge(w,cur_match + 1) == R[j + cur_match]){  //1-based indexing on edges
                    ++cur_match;
                }
                if (cur_match == 0 || cur_match < edge_length){
                    marked[suff_tree.id(v)][0] = 1;
                    marked[suff_tree.id(v)][1] = i;
                    marked[suff_tree.id(v)][2] = cur_match; // Offset on the edge = number of matched chars on it
                    marked[suff_tree.id(v)][3] = v_depth;
                    ++i;
                    j = i + marked[suff_tree.id(v)][3] + cur_match;
                    v = suff_tree.sl(v);
                    break;
                } else{
                    j+= cur_match;
                    v = w;
                }
            }
        }
        std::vector<size_t> table_of_best = calculate_table_of_deepest_nodes(suff_tree, marked);
        int p = 0; 
        while (p + global_offset < block.size()-1){
            size_t idv = table_of_best[p+global_offset];
            size_t length_of_factor = marked[idv][2] + marked[idv][3];
            size_t start_of_factor = marked[idv][1];
            factors.push_back(std::make_tuple(p+global_offset,start_of_factor,length_of_factor));
            p += length_of_factor;
        }
        // Here we assume that have found at least one factor. (Meaning that forall c\in ALPHABET, c is contained in R)
        size_t x = factors.back()[0]; 
        while (x + factors.back()[2]; == block_size){
            b += block_size;
            block_size = std::min(n-b, tau);
            block = T.substr(b, block_size);
            size_t tmp1, tmp2;
            pm_rlz PM{factors.back()[1], factors.back()[2], R, block};
            auto [tmp1, tmp2] = PM.longest_prefix_match();
            factors.back()[1] = tmp1;
            factors.back()[2] = tmp2;
        }
        global_offset = factors.back()[2] - (b - x);
    }

}

// Resulting complexity O(m + tau*log(tau)) = O(m), if tau is chosen properly
