
#pragma once

vector<int64_t> get_border_array(const char* P, int64_t P_len){
    vector<int64_t> B(P_len+1);
    if(B.size() == 1) return B;
    for(int64_t L = 2; L <= P_len; L++){
        int64_t prev_longest = L-1;
        while(prev_longest != 0){
            prev_longest = B[prev_longest];
            if(P[L-1] == P[prev_longest]){
                B[L] = prev_longest + 1;
                break;
            }
        }
    }
   return B;
   
}

// P: search pattern
// P_len: search pattern length
// S: string to search from
// S_len the length of S
// border_array: border array of P
// Returns the number of matches
int64_t kmp(const char* P, int64_t P_len, const char* S, int64_t S_len, const vector<int64_t>& border_array){
    int64_t length = 0;
    int64_t matches = 0;
    for(int64_t i = 0; i < S_len; i++){
        while(length == P_len || P[length] != S[i]){
            length = border_array[length];
            if(length == 0) break;
        }
        if(S[i] == P[length]) length++;
        if(length == P_len) matches++;
    }
    return matches;
}

// TODO: untested for patterns that have length that is not a multiple of 8
int64_t bit_parallel_matching(const char* P, const int64_t P_len, const char* S, const int64_t S_len){
    int64_t n_matches = 0;
    for(int64_t start = 0; start < S_len - P_len + 1; start++){
        bool good = true;
        int64_t n_words = P_len/8;
        for(int64_t w = 0; w < n_words; w++){ // Word-parallel 8 bytes at a time
            good &= (*reinterpret_cast<const uint64_t*>(P + w*8) == 
                     *reinterpret_cast<const uint64_t*>(S + start + w*8));
        }

        // Do the rest
        if(P_len % 8 != 0){
            int64_t rem = P_len - n_words*8; // Number of bytes left over
            good &= *reinterpret_cast<const uint64_t*>(P + n_words*8) >> ((8 - rem)*8) == 
                    *reinterpret_cast<const uint64_t*>(S + start + n_words*8) >> ((8 - rem)*8);
        }
        n_matches += good;
    }
    return n_matches;
}
