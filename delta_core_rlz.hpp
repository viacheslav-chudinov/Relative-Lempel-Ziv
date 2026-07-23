#include "rads_ds.hpp"
#include <sdsl/suffix_trees.hpp> // This may be underlined red in some code redactors
// Depending on your way of integrating sdsl library. In short, in some cases, it 
// should not mean nor break anything

using rlzfactor_t = std::pair<size_t, size_t>;
using factor_t = std::tuple<size_t, size_t, size_t>;
using table_t = std::vector<std::vector<size_t>>;

using node_type = sdsl::cst_sada<>::node_type;
using factors_t = std::vector<std::tuple<size_t, size_t, size_t>>;

inline constexpr size_t SIGMA = 8;
inline char ALPHABET[SIGMA] = {'a','b','c','d','e','f','g', '\0'};

struct Breakup_result{
    table_t marked;
    std::vector<size_t> table_of_best;
};
// So far we implement the case when rads = rads_primitive, as access and updates aare slightly different in the general case
template <typename rads_t>
struct Delta_core_t{
    int l_pow;  // indicating 3*delta*2^{l_pow+1} core-matching query, i.e. the longest
    const std::string_view R;
    const std::string_view T;
    std::vector<rlzfactor_t> factors;
    std::vector<factor_t> pending_factors;
    rads_t main;
    rads_t supplementary;
    std::string batch;
    size_t batch_count;
    size_t batch_capacity = (delta+1)*r << 1;
    std::vector<bool> lookup_new_factors;  // when we factorized eveything till i, new factors maybe in
    Delta_core_t(std::string& ref, std::string& text, size_t tau, size_t delta){
        this->R = ref;
        this->T = text;
        this->l_pow = -1;
        this->main = rads_t rads(R, T, tau);
        this->supplementary = rads_t rads(R, T, tau);
        size_t r = tau / delta + (tau % delta != 0 ? 1 : 0);
        std::string tmp;
        this->batch_capacity = (delta+1)*r << 1;
        tmp.resize(this->batch_capacity, '\0');// Basically it means we're ready to compute local rlz of batches, in which
        // we have 2r phrases of size each of length delta. +1 appears as we'll need toappend to each phrase a terminal symbol to factorize by means of a suffix tree
        this->batch = std::move(tmp);
    }
    bool is_supported_by_main(size_t i){
        if (i < main.buffer_start + main.buffer_size && 
            i >= main.buffer_start - main.l_match){
            return true;
        } else { 
            return false; 
        }
    }
    bool is_supported_by_suppl(size_t i){
        if (i < supplementary.buffer_start + supplementary.buffer_size && 
            i >= supplementary.buffer_start - supplementary.l_match){
            return true;
        } else { 
            return false; 
        }
    }
    char access_T(size_t i){
        if (is_supported_by_main(i)){
            return main.text_at(i - (main.buffer_start - main.l_match));
        } else if (is_supported_by_suppl(i)){
            return supplementary.text_at(i - (supplementary.buffer_start - supplementary.l_match));
        } else{
            throw std::runtime_error("Unsupported index");
        }
    }
    char ref_at(size_t j){
        if (j < R.size()){
            return R[j];
        } else{
            throw std::runtime_error("Bad indexing when accessing the reference ... Access denied");
        } 
    }
    void deposit_to_delta_buffer(size_t b){
        size_t u;
        for (u = 0; u < std::min(delta, n-b); ++u){
            delta_buffer += access_T(b+u);
            delta_buffer += '\0';
        }
        batch_count += u;
        l_pow = -1;
    }
    std::tuple<size_t, size_t, size_t> UpdateMS(size_t l, size_t s, size_t p, size_t query_st = 0);
    std::pair<size_t, size_t> longest_prefix_match(size_t query_st = 0, size_t query_len = 0);
    std::pair<size_t, size_t> core_matching_query(size_t query_st); // answer core-matching query of the form T[query_st...query_st+query_len*3]

    //API from erlz_simple. Does the same thing:
    std::vector<size_t> calculate_table_of_deepest_nodes(sdsl::cst_sada<>& cst, table_t& marked);
    void write_into_marked_1(size_t& id_v, size_t& v_depth, size_t& i, table_t& marked);
    void write_into_marked_0(size_t& id_w, size_t& v_depth, size_t& i, size_t cur_match, table_t& marked);
    void extend(node_type& v, size_t& i, table_t& marked, sdsl::cst_sada<>& suff_tree, size_t& m);

    void decompose_by_metaphrases();
    Breakup_result breakup_delta_buffer();
    void send_factors_out(Breakup_result&& state);

};

template <typename rads_t>
std::tuple<size_t, size_t, size_t> Delta_core_t<rads_t>::UpdateMS(size_t l, size_t s, size_t p, size_t query_st = 0){
    if (l == 0){
        return {1,0,1};
    }
    size_t i = l;
    while (i < l + 1){
        if (access_T(query_st+i-p) < access_T(query_st+i)){               
            i = i - (i>=s ? (i-s)%p : p-((s-i)%p));  
            s = i;                                   
            p = 1;                                   
        } else if(access_T(query_st+i-p) > access_T(query_st+i)){
            p = i - s + 1;
        }
        i++;
    }
    return {l+1, s, p};
}

template <typename rads_t>
std::pair<size_t, size_t> Delta_core_t<rads_t>::longest_prefix_match(size_t query_st = 0, size_t query_len = 0){
    size_t i = 0, p = 0, s = 0, l = 0, l_max = 0;
    size_t m = query_len == 0 ? l_match + buffer_size : query_len;
    size_t n = R.size(); 
    size_t ans_pos = 0;
    while (i < n){
        while (i+l<n && l<m && ref_at(i+l)==access_T(query_st+l)){
            std::tie(l, s, p) = UpdateMS(l, s, p);
        }
        if (l >= l_max){
            ans_pos = i;
            l_max = l;
        }
        bool s_match = true;
        for (size_t b = 0; b < s; ++b){ 
            if(access_T(query_st+b)!=access_T(query_st+p+b)){
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

template <typename rads_t>
std::pair<size_t, size_t> Delta_core_t<rads_t>::core_matching_query(size_t query_st){ // query_st in [0..n-1]
    size_t l, i;
    size_t query_length = (3*delta)<<(l_pow + 1);
    std::tie(i,l) = longest_prefix_match(query_st - (main.buffer_start-main.l_match), query_length);
    if (l != query_length/3){
        return std::make_pair(i, 0);
    } else{
        return std::make_pair(i, l);
    }
}

template <typename rads_t>
std::vector<size_t> Delta_core_t<rads_t>::calculate_table_of_deepest_nodes(sdsl::cst_sada<>& cst, table_t& marked)
{
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
        size_t tmp = marked[id_v][2] + marked[id_v][3];
        if (tmp > val){
            candidate = v;
            val = tmp;
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

template <typename rads_t>
void Delta_core_t<rads_t>::write_into_marked_1(size_t& id_v, size_t& v_depth, size_t& i, table_t& marked)
{   
    marked[id_v][0] = 1;
    marked[id_v][1] = i; // The starting position of the suffix of R concerned rn
    marked[id_v][2] = 0; // Offset on the edge = number of matched chars on it
    marked[id_v][3] = v_depth;
}

template <typename rads_t>
void Delta_core_t<rads_t>::write_into_marked_0(size_t& id_w, size_t& v_depth, size_t& i, size_t cur_match, table_t& marked)
{
    if (marked[id_w][0] == 0 && cur_match >= marked[id_w][2]){
        marked[id_w][1] = i; // The starting position of the suffix of R concerned rn
        marked[id_w][2] = cur_match; // Offset on the edge = number of matched chars on it
        marked[id_w][3] = v_depth; // depth of the anscestor
    }
}


template <typename rads_t>
void Delta_core_t<rads_t>::extend(node_type& v, size_t& i, table_t& marked, sdsl::cst_sada<>& suff_tree, size_t& m)
{
    while (i < m){
        size_t id_v = suff_tree.id(v);
        size_t v_depth = suff_tree.depth(v);
        if (suff_tree.is_leaf(v)){
            v_depth -= 1;
        }
        write_into_marked_1(id_v, v_depth, i, marked);
        i++;
        v = suff_tree.sl(v);
    }
    return ;
}

template <typename rads_t>
void Delta_core_t<rads_t>::breakup_delta_buffer()
{
    size_t m = R.size();
    factors_t factors;
    size_t b = 0;
    size_t block_size;
    block_size = batch_count;
    sdsl::cst_sada<> suff_tree;
    sdsl::construct_im(suff_tree, batch, 1);
    table_t marked(suff_tree.nodes(), std::vector<size_t> (4));
    size_t i = 0;
    size_t j = 0;
    node_type v = suff_tree.root();
    while (i < m){
        while (j < m){
            if (suff_tree.is_leaf(v)){
                size_t id_v = suff_tree.id(v);
                size_t v_depth = suff_tree.depth(v) - 1;
                write_into_marked_1(id_v, v_depth, i, marked);
                v = suff_tree.sl(v);
                i++;
                j = i;
                break;
            }
            // Trying to read a new char
            char c = ref[j];
            node_type w = suff_tree.child(v, c); 
            size_t cur_match = 0;
            size_t v_depth = suff_tree.depth(v);
            int edge_length = suff_tree.depth(w) - v_depth;
            if (edge_length <= 0){
                size_t id_v = suff_tree.id(v);
                write_into_marked_1(id_v, v_depth, i, marked);
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
                size_t id_w = suff_tree.id(w);
                size_t id_v = suff_tree.id(v);
                write_into_marked_0(id_w, v_depth, i, cur_match, marked);
                write_into_marked_1(id_v, v_depth, i, marked);
                j = std::max(i + marked[suff_tree.id(v)][3], i + 1);
                ++i;
                v = suff_tree.sl(v);
                break;
            } else{
                j+= edge_length;
                v = w;
                if (j >= m){
                    extend(v, i, marked, suff_tree, m);
                    break;
                }
            }
        }
    }
    std::vector<size_t> table_of_best = calculate_table_of_deepest_nodes(suff_tree, marked);
    return Breakup_result {std::move(marked), std::move(table_of_best)};
}

// TODO: ADD MERGE SORT LOGIC TO SENDING OUT FACTORS
template <typename rads_t>
void Delta_core_t<rads_t>::send_factors_out(Breakup_result&& state){
    size_t delta_phrase_ind = 0;
    size_t factor_from_cm_ind = 0;
    for (size_t x = 0; x < lookup_new_factors.size(); ++x){
        bool flag = lookup_new_factors[p];  // Becareful, don't switch the type to auto. Vector of bools is meemory-packed and returns a proxy-type
        if (flag){
            size_t p = delta_phrase_ind; 
            size_t idv = state.table_of_best[p];
            size_t length_of_factor = state.marked[idv][2] + state.marked[idv][3];
            size_t start_of_factor = state.marked[idv][1];
            factors.push_back(std::make_tuple(p,start_of_factor, length_of_factor));
            p += length_of_factor;
            while (ref_at(p-1)!='\0' && p < batch_count){
                idv = state.table_of_best[p];
                length_of_factor = state.marked[idv][2] + state.marked[idv][3];
                start_of_factor = state.marked[idv][1];
                factors.push_back(std::make_tuple(p,start_of_factor, length_of_factor));
                p += length_of_factor;
            }
            factors.back[2] -= 1;
            delta_phrase_ind = p + 1;
        } else{
            factors.push_back(pending_factors[factor_from_cm_ind]);
            factor_from_cm_ind++;
        }
    }
    size_t p = 0; 
    while (p < block_size){
        size_t idv = state.table_of_best[p];
        size_t length_of_factor = state.marked[idv][2] + state.marked[idv][3];
        size_t start_of_factor = state.marked[idv][1];
        factors.push_back(std::make_tuple(p,start_of_factor, length_of_factor));
        p += length_of_factor;
    }
}

template <typename rads_t>
void Delta_core_t<rads_t>::decompose_by_metaphrases(){
    size_t n = T.size();
    size_t b = 0;
    size_t i, l;
    main.fill();
    supplementary.fill();
    while(b < n){
        bool can_extend = true;
        while(can_extend){
            size_t qr_step = delta<<(l_pow+1);
            bool can_run_next_cm_query = is_supported_by_main(b + (3*qr_step)-1);
            bool lost_access = false;
            lookup_new_factors.push_back(0);
            while ((!can_run_next_cm_query) && (!lost_access)){
                main.update();
                can_run_next_cm_query = is_supported_by_main(b + (3*qr_step)-1);
                if (is_supported_by_main(b)){
                    supplementary.update();
                } else{
                    lost_access = true;
                    break;
                }
            }
            if (!is_supported_by_main(b) && !is_supported_by_suppl(b)) {
                throw std::runtime_error("Lost access to meta-phrase start ... Execution failed");
            }
            if (lost_access){
                deposit_to_delta_buffer(b);
                can_extend = false;
                break;
            }
            std::tie(i,l) = core_matching_query(b);
            if (l == 0){
                deposit_to_delta_buffer(b);
                can_extend = false;
                break;
            }
            lookup_new_factors.push_back(1);
            l_pow++;
            pending_factors.push_back(std::make_tuple(b+qr_step, i, l)); // PONDER A CHANGE TO RLZFACTOR_T 
        }
        b += delta<<(l_pow+1);
        if (batch_count >= batch_capacity){
            Breakup_result state = breakup_delta_buffer();
            send_factors_out(std::move(state));
        }
    }
    if (batch_count > 0){
        Breakup_result state = breakup_delta_buffer();
        send_factors_out(std::move(state));
    } else if (pending_factors.size() > 0){
        for (auto s : pending_factors){
            factors.push_back(s);
        }
    }
}

