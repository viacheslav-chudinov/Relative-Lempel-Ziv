#include "rads.hpp"
#include <sdsl/suffix_trees.hpp> // This may be underlined red in some code redactors
// Depending on your way of integrating sdsl library. In short, in some cases, it 
// should not mean nor break anything

using factor_t = std::pair<size_t, size_t>;
using table_t = std::vector<std::vector<size_t>>;

using node_type = sdsl::cst_sada<>::node_type;
using factors_t = std::vector<std::tuple<size_t, size_t, size_t>>;

inline constexpr size_t SIGMA = 8;
inline char ALPHABET[SIGMA] = {'a','b','c','d','e','f','#', '\0'};

struct Breakup_result{
    table_t marked;
    std::vector<size_t> table_of_best;
};
// So far we implement the case when rads = rads_primitive, as access and updates are slightly different in the general case
template <typename rads_t>
struct Delta_core_t{
    int l_pow;  // indicating 3*delta*2^{l_pow+1} core-matching query, i.e. the longest
    const std::string_view R;
    const std::string_view T;
    std::vector<factor_t> factors;
    std::vector<factor_t> pending_factors;
    rads_t main;
    rads_t supplementary;
    std::string batch;
    size_t delta;
    size_t batch_count;
    size_t batch_capacity;
    size_t main_steps;
    size_t suppl_steps;
    std::vector<bool> lookup_new_factors;  // when we factorized eveything till i, new factors maybe in
    Delta_core_t(std::string& ref, std::string& text, size_t tau, size_t delta) : l_pow(-1),  R(ref), T(text),
    main(ref, text, tau), supplementary(ref, text, tau), delta(delta), batch_count(0), main_steps(0), suppl_steps(0){
        size_t r = tau / delta + (tau % delta != 0 ? 1 : 0);
        std::string tmp;
        this->batch_capacity = (delta+1)*r << 1;
        tmp.resize(this->batch_capacity, '#');// Basically it means we're ready to compute local rlz of batches, in which
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
            std::cout << "Unsupported index i = " << i << std::endl;
            throw std::runtime_error("Unsupported index i");
        }
    }
    char ref_at(size_t j){
        if (j < R.size()){
            return R[j];
        } else{
            throw std::runtime_error("Bad indexing when accessing the reference ... Access denied");
        } 
    }
    size_t deposit_to_delta_buffer(size_t b, size_t u = 0){
        // std::cout << "We are in deposit_to_delta_buffer" << std::endl;
        size_t n = T.size();
        while (u < std::min(delta, n-b) && (is_supported_by_main(b+u) || is_supported_by_suppl(b+u))){
            batch[batch_count+u] = access_T(b+u);
            u++;
        }
        return u;
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
std::tuple<size_t, size_t, size_t> Delta_core_t<rads_t>::UpdateMS(size_t l, size_t s, size_t p, size_t query_st){
    // std::cout << "We are in UpdateMS" <<std::endl;
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
std::pair<size_t, size_t> Delta_core_t<rads_t>::longest_prefix_match(size_t query_st, size_t query_len){
    // std::cout << "We are in longest_prefix_match" <<std::endl;
    size_t i = 0, p = 0, s = 0, l = 0, l_max = 0;
    size_t m = query_len == 0 ? main.l_match + main.buffer_size : query_len;
    size_t n = R.size(); 
    size_t ans_pos = 0;
    while (i < n){
        while (i+l<n && l<m && ref_at(i+l)==access_T(query_st+l)){
            std::tie(l, s, p) = UpdateMS(l, s, p, query_st);
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
    // std::cout << "We are in core_matching_query with query_st = " << query_st <<std::endl;
    size_t l, i;
    size_t query_length = delta<<(l_pow + 1);
    std::tie(i,l) = longest_prefix_match(query_st + query_length, query_length);
    if (l != query_length){
        return std::make_pair(i, 0);
    } else{
        return std::make_pair(i, l);
    }
}

template <typename rads_t>
std::vector<size_t> Delta_core_t<rads_t>::calculate_table_of_deepest_nodes(sdsl::cst_sada<>& cst, table_t& marked)
{
    // std::cout << "We are in calculate_table_of_best" <<std::endl;
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
Breakup_result Delta_core_t<rads_t>::breakup_delta_buffer()
{
    // std::cout << "We are in breakup_delta_buffer" << std::endl;
    size_t m = R.size();
    size_t b = 0;
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

            size_t id_v = suff_tree.id(v);
            size_t v_depth  = suff_tree.depth(v);
            write_into_marked_1(id_v, v_depth, i, marked);

            // Trying to read a new char
            char c = ref_at(j);
            node_type w = suff_tree.child(v, c); 
            size_t cur_match = 0;
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
            while (cur_match < edge_length && j + cur_match < m && suff_tree.edge(w,v_depth + cur_match + 1) == ref_at(j + cur_match)){  //1-based indexing on edges
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
                j += edge_length;
                size_t id_w = suff_tree.id(w);
                size_t w_depth = v_depth + edge_length;    // == suff_tree.depth(w)
                write_into_marked_1(id_w, w_depth, i, marked);
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

template <typename rads_t>
void Delta_core_t<rads_t>::send_factors_out(Breakup_result&& state){
    // std::cout << "We are in send_factors_out" << std::endl;
    size_t delta_phrase_ind = 0;
    size_t factor_from_cm_ind = 0;
    // std::cout << "lookup_new_factors trace is " << std::endl;
    // for (int i = 0; i < lookup_new_factors.size(); ++i){
    //     std::cout << lookup_new_factors[i] << std::endl;
    // } 
    // std::cout << "pending factors " << std::endl;
    // for (auto s : pending_factors){
    //     std::cout << s.first << " " << s.second << std::endl;
    // }
    // std::cout << "batch_count = " << batch_count << std::endl;
    // std::cout << "batch is : \n" << batch << std::endl;
    for (size_t x = 0; x < lookup_new_factors.size(); ++x){
        bool flag = lookup_new_factors[x];  // Becareful, don't switch the type to auto. Vector of bools is memory-packed and returns a proxy-type
        // std::cout << "flag = " << flag <<std::endl;
        if (!flag){
            size_t p = delta_phrase_ind; 
            size_t idv = state.table_of_best[p];
            size_t length_of_factor = state.marked[idv][2] + state.marked[idv][3];
            size_t start_of_factor = state.marked[idv][1];
            factors.push_back(std::make_pair(start_of_factor, length_of_factor));
            // std::cout << "p = " << p << std::endl;
            // std::cout << "First factor is " << factors.back().first << " " << factors.back().second <<std::endl;
            // std::cout << "First factor is " << R[factors.back().first] <<std::endl;
            p += length_of_factor;
            while (p < batch_count && batch[p]!='#'){  // MAY BE SIMPLIFIED TO batch[p-1]!='#'
                // std::cout << "batch[p] = " << batch[p] <<std::endl;
                idv = state.table_of_best[p];
                // std::cout << "p=" << p << " idv=" << idv
                // << " off=" << state.marked[idv][2]
                // << " dep=" << state.marked[idv][3]
                // << " start=" << state.marked[idv][1] << std::endl;
                length_of_factor = state.marked[idv][2] + state.marked[idv][3];
                if (length_of_factor == 0){
                    throw std::runtime_error("Zero-length factor at batch position " + std::to_string(p));
                }
                start_of_factor = state.marked[idv][1];
                factors.push_back(std::make_pair(start_of_factor, length_of_factor));
                p += length_of_factor;
            }
            delta_phrase_ind = p + 1;
        } else{
            // std::cout << "Adding a cm factor " <<std::endl;
            factors.push_back(pending_factors[factor_from_cm_ind]);
            factor_from_cm_ind++;
        }
    }
    lookup_new_factors.clear();
    pending_factors.clear();
}

template <typename rads_t>
void Delta_core_t<rads_t>::decompose_by_metaphrases(){
    size_t n = T.size();
    size_t b = 0;
    size_t i, l;
    main.fill();
    supplementary.fill();
    while (b < n){
        lookup_new_factors.push_back(0);
        l_pow = -1;
        if (supplementary.buffer_start > b + supplementary.l_match){
            throw std::runtime_error("Supplementary overshot meta-phrase start b");
        }
        if ((main.buffer_start + main.buffer_size - 1) < b){
            while (!is_supported_by_main(b)){
                if (main.buffer_start >= n){
                    throw std::runtime_error("Main exhausted while catching up to b");
                }
                main.update();
            }
        }
        while (!is_supported_by_suppl(b)){
            if (supplementary.buffer_start >= n){
                throw std::runtime_error("Supplementary exhausted before reaching b");
            }
            supplementary.update();
        }
        bool can_extend = true;
        while (can_extend){
            size_t qr_step = delta << (l_pow + 1);
            if (b + 3*qr_step > n){
                can_extend = false;
                break;
            }
            bool lost_access = !is_supported_by_main(b);
            while (!is_supported_by_main(b + 3*qr_step - 1) && !lost_access){
                if (main.buffer_start >= n){
                    lost_access = true;
                    break;
                }
                main.update();
                if (!is_supported_by_main(b)){
                    lost_access = true;
                    break;
                }
                supplementary.update();
            }
            if (lost_access || !is_supported_by_main(b)){
                can_extend = false;
                break;
            }
            std::tie(i, l) = core_matching_query(b);
            if (l == 0){
                can_extend = false;
                break;
            }
            lookup_new_factors.push_back(1);
            l_pow++;
            pending_factors.push_back(std::make_pair(i, l));
        }
        size_t need = std::min(delta, n - b);
        size_t u = deposit_to_delta_buffer(b, 0);
        while (u < need){
            if (supplementary.buffer_start >= n){
                throw std::runtime_error("Supplementary exhausted mid delta-prefix deposit");
            }
            supplementary.update();
            u = deposit_to_delta_buffer(b, u);
        }
        if (u != need){
            throw std::runtime_error("Incomplete delta-prefix deposit");
        }
        batch_count += u;
        batch[batch_count] = '#';
        batch_count += 1;

        b += delta << (l_pow + 1);

        if (batch_count >= batch_capacity){
            Breakup_result state = breakup_delta_buffer();
            send_factors_out(std::move(state));
            batch_count = 0;
        }
    }

    if (batch_count > 0){
        Breakup_result state = breakup_delta_buffer();
        send_factors_out(std::move(state));
        batch_count = 0;
    } else if (pending_factors.size() > 0){
        for (const auto& s : pending_factors){
            factors.push_back(s);
        }
        pending_factors.clear();
    }
}