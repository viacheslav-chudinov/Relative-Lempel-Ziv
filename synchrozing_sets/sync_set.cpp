#include "sync_set.hpp"
#include "../compute_period.hpp"
#include "../hashes/karp-rabin.hpp"
#include <algorithm>

uint32_t eval_hash(uint128_t x, std::vector<uint32_t>& coeffs, uint32_t d_){
    using boost::multiprecision::uint256_t;
    uint256_t res = static_cast<uint256_t>(0);
    uint256_t P = static_cast<uint256_t>(d_);
    for (long int i = coeffs.size() - 1; i >= 0; --i)
        res = (((res * x)%P) + coeffs[i])%P;
    return static_cast<uint128_t>(res);
}

bool Sync_Set::sync(const uint32_t ind, const uint128_t kr_fingerprint, const access_func& access){
    uint32_t per_plus = compute_period(ind + 1, tau_<<1, access);
    uint32_t per_minus = compute_period(ind, (tau_<<1) - 1, access);
    if (per_plus > tau_/3 && per_minus > tau_/3){
        return eval_hash(kr_fingerprint, hash_coeffs, d_) == 0;
    }
    uint32_t per = compute_period(ind, tau_<<1, access);
    if (per > tau_/3){
        return true;
    } 
    return false;
}

void Sync_Set::enum_block(const uint32_t ind, const access_func& access, std::vector<uint32_t>& sync_pos){
    uint32_t len = 3*tau_ - 1;
    uint32_t period = compute_period(ind + tau_ - 1, tau_ + 1, access); // covers [tau-1, 2tau-1]
    long int l = this->tau_;
    uint32_t r = this->tau_;
    if (period <= tau_/3){
        long int l_prime = (long int)tau_ - 2;
        uint32_t r_prime = (tau_<<1);
        for (l_prime = (long int)tau_ - 2; l_prime >= 0; l_prime--){
            if (access(ind + l_prime) != access(ind + l_prime + period)){
                l_prime++;
                break;
            }
        }
        if (l_prime < 0) l_prime = 0;
        for (r_prime = (tau_<<1); r_prime < len; r_prime++){
            if (access(ind + r_prime) != access(ind + r_prime - period)){
                break;
            }
        }
        if (((long)r_prime - l_prime)>=(2*tau_ - 1)){
            l = l_prime - 1;
            r = r_prime - (tau_<<1) + 1;
        }
    }
    uint128_t base_powed = compute_pow_modulo_n_long(this->base_, (tau_<<1) - 1, this->q_prime);
    uint128_t kr_fp = 0;
    if (l>0){
        kr_fp = phi_at_pos(ind, tau_<<1, this->base_, this->q_prime, access);
        if (eval_hash(kr_fp, this->hash_coeffs, this->d_) == 0){
            sync_pos.push_back(ind);
        }
    }
    for (int s = 1; s < l; s++){
        kr_fp = phi_step(ind + s, tau_<<1, this->base_, base_powed, kr_fp, this->q_prime, access);
        if (eval_hash(kr_fp, this->hash_coeffs, this->d_) == 0){
            sync_pos.push_back(static_cast<uint32_t>(ind+s));
        }
    }
    if (l >= 0 && l < this->tau_){
        sync_pos.push_back(ind+l);
    }
    if (r < this->tau_){
        sync_pos.push_back(ind+r);
    }
    if (r < (this->tau_ - 1)){
        kr_fp = phi_at_pos(ind + r + 1, tau_<<1, base_, q_prime, access);
        if (eval_hash(kr_fp, this->hash_coeffs, this->d_) == 0){
            sync_pos.push_back(static_cast<uint32_t>(ind+r+1));
        }
    }
    for (uint32_t s = r + 2; s < this->tau_; s++){
        kr_fp = phi_step(ind + s, tau_<<1, this->base_, base_powed, kr_fp, this->q_prime, access);
        if (eval_hash(kr_fp, this->hash_coeffs, this->d_) == 0){
            sync_pos.push_back(ind+s);
        }
    }
}

std::vector<uint32_t> Sync_Set::enumerate_sync(const uint32_t ind, const uint32_t length, const access_func& access){
    if (length < (tau_<<1)){
        throw std::runtime_error("Too short to admit a sync position");
    } else if (length < (3*tau_-1)){
        std::vector<uint32_t> sync_pos;
        access_func ultimate_access = [ind, length, &access](uint32_t jj) -> char {
            if ((jj - ind) < length){
                return access(jj);
            }
            return '#';
        };
        enum_block(ind, ultimate_access, sync_pos);
        auto it = std::lower_bound(sync_pos.begin(), sync_pos.end(), ind + length - (tau_<<1) + 1);
        size_t idx = it - sync_pos.begin();
        sync_pos.resize(idx);
        return sync_pos;
    } else{
        std::vector<uint32_t> sync_pos;
        uint32_t h = 0;
        uint32_t right = (length-3*tau_+1)/tau_;
        for (h = 0; h <= right; ++h){
            enum_block(ind + (h*tau_), access, sync_pos);
        }
        std::vector<uint32_t> sync_pos1;
        enum_block(ind + (length - (3*tau_) + 1), access, sync_pos1);
        if (!sync_pos1.empty()){
            uint32_t xx = 0;
            if (sync_pos.empty()){
                return sync_pos1;
            }
            uint32_t target = sync_pos.back();
            while (xx < sync_pos1.size() && sync_pos1[xx] <= target){
                xx++;
            }
            while (xx < sync_pos1.size()){
                sync_pos.push_back(sync_pos1[xx]);
                xx++;
            }
        }
        return sync_pos;
    }
}


bool Sync_Set::is_sparse_and_dense(const uint32_t ind, const uint32_t len, const access_func& access){
    if (len < 3*tau_ - 1){
        return true;
    }
    uint32_t sparsity_bound = 72*k_+2;
    std::vector<uint32_t> sync_positions = enumerate_sync(ind, len, access);
    const uint32_t last_i = ind + len - 3*tau_ + 1;
    uint32_t left = 0, right = 0;
    for (uint32_t i = ind; i <= last_i; ++i){
        while (right < sync_positions.size() && sync_positions[right] < i + tau_){
            ++right;
        }
        while (left < right && sync_positions[left] < i){
            ++left;
        }
        const uint32_t count = right - left;
        if (count > sparsity_bound){
            return false;
        }
        if (count == 0){
            uint32_t per = compute_period(i, 3*tau_ - 1, access);
            if (per > tau_/3){
                return false;
            }
        }
    }
    return true;
}

void Sync_Set::reinit(){
    this->base_ = random_in_range_prime(this->q_prime);
    if (tau_ >= 12*k_){
        for (uint32_t i = 0; i < k_; ++i){
            this->hash_coeffs[i] = random_in_range_prime_with_0(d_);
        }
    }
}