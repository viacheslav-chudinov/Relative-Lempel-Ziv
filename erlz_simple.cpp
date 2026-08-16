#include <stack>
#include "pm&rlz_utils.hpp"

// TODO :
// Add std::move logic, because currently std::strings 
// are copied each time PM is invoked

using table_t = std::vector<std::vector<uint_t>>;

using node_type = sdsl::cst_sada<>::node_type;



std::vector<uint_t> erlz_simple::calculate_table_of_deepest_nodes(sdsl::cst_sada<>& cst, table_t& marked)
{   
    // Post-order traversal
    node_type root = cst.root();
    std::stack<node_type> stack1;
    std::vector<std::pair<node_type, uint_t>> stack2;
    stack1.push(cst.root());
    stack2.push_back(std::make_pair(cst.root(), NONE));
    while(!stack1.empty()){
        node_type v = stack1.top();
        stack1.pop();
        for (const auto& w : cst.children(v)){
            stack2.push_back(std::make_pair(w, NONE));
            stack1.push(w);
        }
    }

    long int g = stack2.size() - 1;
    while(g >= 0){
        node_type v = stack2[g].first;
        uint_t id_v = cst.id(v);
        uint_t suff_ind_of_R = NONE;
        if (marked[id_v][0]==0 && marked[id_v][2] == 0){
            for (const auto& w : cst.children(v)){
                if(marked[cst.id(w)][1] < NONE) suff_ind_of_R = marked[cst.id(w)][1];
            }
            if (suff_ind_of_R != NONE){
                marked[id_v][1] = suff_ind_of_R;
                marked[id_v][3] = cst.depth(v);
            }
        }
        g--;
    }
    stack2.clear();

    // Pre-order traversal:

    std::vector<uint_t> table_of_best(cst.csa.size());
    std::stack<std::tuple<node_type, uint_t, node_type>> stack;
    stack.push({root, 0, root});
    while(!stack.empty()){
        std::tuple state = stack.top();
        node_type v = std::get<0>(state);
        uint_t val = std::get<1>(state);
        node_type candidate = std::get<2>(state);
        stack.pop();
        uint_t id_v = cst.id(v);
        uint_t tmp = marked[id_v][2] + marked[id_v][3];
        if (tmp > val){
            candidate = v;
            val = tmp;
        }
        if (cst.is_leaf(v)){
            table_of_best[cst.sn(v)] = cst.id(candidate);
        } else {
            for (const auto& w : cst.children(v)){
                stack.push({w, val, candidate});
            }
        }
    }
    return table_of_best;
}


void erlz_simple::write_into_marked_1(uint_t& id_v, uint_t& v_depth, uint_t& i, table_t& marked){
    marked[id_v][0] = 1;
    marked[id_v][1] = i; // The starting position of the suffix of R concerned rn
    marked[id_v][2] = 0; // Offset on the edge = number of matched chars on it
    marked[id_v][3] = v_depth;
}
void erlz_simple::write_into_marked_0(uint_t& id_w, uint_t& v_depth, uint_t& i, uint_t cur_match, table_t& marked){
    if (marked[id_w][0] == 0 && cur_match >= marked[id_w][2]){
        marked[id_w][1] = i; // The starting position of the suffix of R concerned rn
        marked[id_w][2] = cur_match; // Offset on the edge = number of matched chars on it
        marked[id_w][3] = v_depth; // depth of the anscestor
    }
}

void erlz_simple::extend(node_type& v, uint_t& i, table_t& marked, sdsl::cst_sada<>& suff_tree, uint_t& m){
    while (i < m){
        uint_t id_v = suff_tree.id(v);
        uint_t v_depth = suff_tree.depth(v);
        if (suff_tree.is_leaf(v)){
            v_depth -= 1;
        }
        write_into_marked_1(id_v, v_depth, i, marked);
        i++;
        v = suff_tree.sl(v);
    }
    return ;
}


factors_t erlz_simple::rlz_compress()
{
    uint_t m = ref.size();
    uint_t n = text.size();
    uint_t global_offset = 0;
    factors_t factors;
    uint_t b = 0;
    uint_t block_size;
    while (b < n){
        block_size = std::min(tau,n-b);
        std::string block = text.substr(b, block_size);
        sdsl::cst_sada<> suff_tree;
        sdsl::construct_im(suff_tree, block, 1);
        uint_t root_id = suff_tree.id(suff_tree.root());
        std::vector<uint_t> init_row = {0, NONE, 0, 0}; 
        table_t marked(suff_tree.nodes(), init_row);
        uint_t i = 0;
        uint_t j = 0;
        node_type v = suff_tree.root();
        while (i < m){
            while (j < m){
                if (suff_tree.is_leaf(v)){
                    uint_t id_v = suff_tree.id(v);
                    uint_t v_depth = suff_tree.depth(v) - 1;
                    write_into_marked_1(id_v, v_depth, i, marked);
                    v = suff_tree.sl(v);
                    j = i + v_depth;
                    i++;
                    break;
                }
                uint_t id_v = suff_tree.id(v);
                uint_t v_depth  = suff_tree.depth(v);
                write_into_marked_1(id_v, v_depth, i, marked);
                // Trying to read a new char
                char c = ref[j];
                node_type w = suff_tree.child(v, c); 
                uint_t cur_match = 0;
                int edge_length = suff_tree.depth(w) - v_depth;
                if (edge_length <= 0){
                    j = std::max(i + marked[suff_tree.id(v)][3], i + 1);
                    ++i;
                    v = suff_tree.sl(v);
                    if (j >= m) 
                        extend(v, i, marked, suff_tree, m);
                    break;
                }
                while (cur_match < edge_length && j + cur_match < m && suff_tree.edge(w,v_depth + cur_match + 1) == ref[j + cur_match]){  //1-based indexing on edges
                    ++cur_match;
                }
                if (cur_match < edge_length){
                    uint_t id_w = suff_tree.id(w);
                    write_into_marked_0(id_w, v_depth, i, cur_match, marked);
                    j = std::max(i + marked[suff_tree.id(v)][3], i + 1);
                    ++i;
                    v = suff_tree.sl(v);
                    break;
                } else{
                    j += edge_length;
                    uint_t id_w = suff_tree.id(w);
                    uint_t w_depth = v_depth + edge_length;    // == suff_tree.depth(w)
                    write_into_marked_1(id_w, w_depth, i, marked);
                    v = w;
                    if (j >= m){
                        extend(v, i, marked, suff_tree, m);
                        break;
                    }
                }
            }
        }
        std::vector<uint_t> table_of_best = calculate_table_of_deepest_nodes(suff_tree, marked);
        int p = 0; 
        while (p + global_offset < block.size()){
            std::cout << "p = " << p << std::endl;
            uint_t idv = table_of_best[p+global_offset];
            uint_t length_of_factor = marked[idv][2] + marked[idv][3];
            uint_t start_of_factor = marked[idv][1];
            factors.push_back(std::make_tuple(p+global_offset,start_of_factor, length_of_factor));
            p += length_of_factor;
        }
        // Here we assume that have found at least one factor. (Meaning that forall c\in ALPHABET, c is contained in R)
        uint_t x = b + std::get<0>(factors.back());
        if (b + block_size >= n) { b = n; }
        while (b < n && b + block_size < n && x + std::get<2>(factors.back()) == (b+block_size)){
            b += block_size;
            block_size = std::min(n-b, tau);
            block = text.substr(b, block_size);
            pm_rlz PM{std::get<1>(factors.back()), std::get<2>(factors.back()), ref, block};
            auto pair = PM.longest_prefix_match();
            std::get<1>(factors.back()) = pair.first;
            std::get<2>(factors.back()) = pair.second;
        }
        global_offset = x + std::get<2>(factors.back()) - b;
    }
    return factors;
}

// Resulting complexity O(m + tau*log(tau)) = O(m), if tau is chosen properly

