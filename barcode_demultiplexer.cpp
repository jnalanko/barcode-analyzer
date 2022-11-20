#include <iostream>
#include <algorithm>
#include "SeqIO.hh"
#include "cxxopts.hpp"


// Table mapping ascii values of characters to their reverse complements,
// lower-case to lower case, upper-case to upper-case. Non-ACGT characters
// are mapped to themselves.
static constexpr unsigned char rc_table[256] =
{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19,
20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37,
38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55,
56, 57, 58, 59, 60, 61, 62, 63, 64, 84, 66, 71, 68, 69, 70, 67, 72, 73,
74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 65, 85, 86, 87, 88, 89, 90, 91,
92, 93, 94, 95, 96, 116, 98, 103, 100, 101, 102, 99, 104, 105, 106, 107,
108, 109, 110, 111, 112, 113, 114, 115, 97, 117, 118, 119, 120, 121, 122,
123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136, 137,
138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149, 150, 151, 152,
153, 154, 155, 156, 157, 158, 159, 160, 161, 162, 163, 164, 165, 166, 167,
168, 169, 170, 171, 172, 173, 174, 175, 176, 177, 178, 179, 180, 181, 182,
183, 184, 185, 186, 187, 188, 189, 190, 191, 192, 193, 194, 195, 196, 197,
198, 199, 200, 201, 202, 203, 204, 205, 206, 207, 208, 209, 210, 211, 212,
213, 214, 215, 216, 217, 218, 219, 220, 221, 222, 223, 224, 225, 226, 227,
228, 229, 230, 231, 232, 233, 234, 235, 236, 237, 238, 239, 240, 241, 242,
243, 244, 245, 246, 247, 248, 249, 250, 251, 252, 253, 254, 255};

static constexpr char get_rc(char c){
    return rc_table[(unsigned char)c];
}

void rc_c_string(char* S, int64_t len){
    std::reverse(S, S + len);
    for(int64_t i = 0; i < len; i++){
        S[i] = get_rc(S[i]);
    }
}

vector<string> read_lines(string filename){
    vector<string> lines;
    throwing_ifstream in(filename);
    string line;
    while(in.getline(line)){
        lines.push_back(line);
    }
    return lines;
}


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

int main(int argc, char** argv){

    cxxopts::Options opts(argv[0], "Barcode demultiplexing");

    opts.add_options()
        ("i", "The sequence file in fasta or fastq format.", cxxopts::value<string>())
        ("b", "A file containing the barcodes, one per line.", cxxopts::value<string>())
        ("v,verbose", "Verbose output.", cxxopts::value<bool>()->default_value("false"))
        ("h,help", "Print usage")
    ;

    auto opts_parsed = opts.parse(argc, argv);

    if (argc == 1 || opts_parsed.count("help")){
        std::cerr << opts.help() << std::endl;
        return 1;
    }

    string seq_file = opts_parsed["i"].as<string>();
    string barcode_file = opts_parsed["b"].as<string>();
    bool verbose = opts_parsed["v"].as<bool>();

    vector<string> barcodes = read_lines(barcode_file);
    vector<string> barcodes_rc = barcodes;
    for(int64_t i = 0; i < barcodes_rc.size(); i++) 
        rc_c_string(barcodes_rc[i].data(), barcodes_rc[i].size());

    vector<vector<int64_t> > barcode_border_arrays;
    for(const string& P : barcodes) barcode_border_arrays.push_back(get_border_array(P.c_str(), P.size()));

    vector<vector<int64_t> > barcode_rc_border_arrays;
    for(const string& P : barcodes_rc) barcode_rc_border_arrays.push_back(get_border_array(P.c_str(), P.size()));

    SeqIO::Reader<> in(seq_file);

    vector<int64_t> barcode_counts(barcodes.size());

    vector<int64_t> barcodes_found; // Used only in verbose mode
    while(true){
        int64_t len = in.get_next_read_to_buffer();
        if(len == 0) break;
        char* seq = in.read_buf;
        
        if(verbose) barcodes_found.clear();    
        for(int64_t barcode_idx = 0; barcode_idx < barcodes.size(); barcode_idx++){
            int64_t occurrences = 0;
            //occurrences += kmp(barcodes[barcode_idx].c_str(), barcodes[barcode_idx].size(), seq, len, barcode_border_arrays[barcode_idx]);
            //occurrences += kmp(barcodes_rc[barcode_idx].c_str(), barcodes_rc[barcode_idx].size(), seq, len, barcode_rc_border_arrays[barcode_idx]);
            occurrences += bit_parallel_matching(barcodes[barcode_idx].c_str(), barcodes[barcode_idx].size(), seq, len);
            occurrences += bit_parallel_matching(barcodes_rc[barcode_idx].c_str(), barcodes_rc[barcode_idx].size(), seq, len);
            barcode_counts[barcode_idx] += occurrences;
            if(verbose && occurrences > 0) barcodes_found.push_back(barcode_idx);
        }

        if(verbose && barcodes_found.size() > 1){
            cout << "Multiple barcodes in a sequence: " << seq << endl;
            cout << "^ Had barcodes:";
            for(int64_t b : barcodes_found) cout << " " << b;
            cout << endl;
        }
    }

    for(int64_t i = 0; i < barcode_counts.size(); i++){
        cout << "Barcode " << i << ": " << barcode_counts[i] << endl;
    }

}