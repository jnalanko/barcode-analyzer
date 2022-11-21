#include <iostream>
#include <algorithm>
#include "SeqIO.hh"
#include "cxxopts.hpp"
#include "aho_corasick.hh"

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

void analyze(const string& seq_file, const string& barcode_file, ostream& output, bool verbose){

    vector<string> barcodes = read_lines(barcode_file);
    int64_t n_barcodes = barcodes.size();

    // Add reverse complements
    for(int64_t i = 0; i < n_barcodes; i++){
        string S = barcodes[i];
        rc_c_string(S.data(), S.size());
        barcodes.push_back(S);
    }
        
    // Build the Aho-Corasick trie
    aho_corasick::trie trie;
    for(const string& B : barcodes) trie.insert(B);

    SeqIO::Reader<> in(seq_file);

    vector<int64_t> global_counts(n_barcodes); // Counts of barcodes in all sequences
    vector<int64_t> local_counts(n_barcodes); // Counts of barcodes in the current sequence

    // List of distinct barcodes found in the current sequence
    vector<int64_t> local_barcodes_found;

    int64_t n_seqs_with_multiple_barcodes = 0;

    while(true){
        int64_t len = in.get_next_read_to_buffer();
        if(len == 0) break;
        char* seq = in.read_buf;
        
        auto AC_result = trie.parse_text(seq);
        for(auto x : AC_result){
            // The modulo is to map the reverse complement barcodes to the same barcode as the original
            int64_t barcode_idx = x.get_index() % n_barcodes;

            //global_counts[barcode_idx]++;
            if(local_counts[barcode_idx] == 0){
                local_barcodes_found.push_back(barcode_idx);
            }
            local_counts[barcode_idx]++;
        }

        if(local_barcodes_found.size() >= 2){
            // Multiple distinct barcodes in this sequence
            n_seqs_with_multiple_barcodes++;
            if(verbose){
                output << "Mixed barcodes in sequence: " << in.header_buf << "\n";
                output << "Found barcodes: ";
                for(int64_t i = 0; i < local_barcodes_found.size(); i++)
                    output << (i == 0 ? "" : " ") << local_barcodes_found[i];
                output << "\n";
            }
        } else{
            // Add local counts to global counts
            for(int64_t x : local_barcodes_found){
                global_counts[x] += local_counts[x];
            }
        }

        // Clear local counters
        for(int64_t x : local_barcodes_found) local_counts[x] = 0;
        local_barcodes_found.clear();
    }

    for(int64_t i = 0; i < global_counts.size(); i++){
        // Print barcodes in 1-based indexing
        output << "Barcode " << i+1 << ": " << global_counts[i] << endl;
    }
    output << "Mixed: " << n_seqs_with_multiple_barcodes << endl;
}

int analyze_main(int argc, char** argv){
    cxxopts::Options opts(argv[0], "Search for barcode sequences inside a fasta/fastq file.");

    opts.add_options()
        ("i", "The sequence file in fasta or fastq format.", cxxopts::value<string>())
        ("o", "Output file. If not given, prints to stdout.", cxxopts::value<string>())
        ("b", "A file containing the barcodes, one per line.", cxxopts::value<string>())
        ("v,verbose", "Verbose output.", cxxopts::value<bool>()->default_value("false"))
        ("h,help", "Print usage")
    ;

    auto opts_parsed = opts.parse(argc, argv);

    if (argc == 1 || opts_parsed.count("help")){
        std::cerr << opts.help() << std::endl;
        return 1;
    }

    bool to_stdout = false;
    string seq_file = opts_parsed["i"].as<string>();
    string barcode_file = opts_parsed["b"].as<string>();
    string output_file;
    try{
        output_file = opts_parsed["o"].as<string>();
    } catch(cxxopts::exceptions::option_has_no_value& e){
        to_stdout = true;
    }
    bool verbose = opts_parsed["v"].as<bool>();

    if(to_stdout) analyze(seq_file, barcode_file, cout, verbose);
    else{
        ofstream out(output_file);
        analyze(seq_file, barcode_file, out, verbose);
    }

    return 0;

}


void filter_barcodes(const string& seq_file, const string& barcode_file, ostream& output){
    
}

int filter_main(int argc, char** argv){
    cxxopts::Options opts(argv[0], "Remove all sequence from a fasta/fastq file that have a barcode sequence.");

    opts.add_options()
        ("i", "The sequence file in fasta or fastq format.", cxxopts::value<string>())
        ("o", "Output file. If not given, prints to stdout.", cxxopts::value<string>())
        ("b", "A file containing the barcodes, one per line.", cxxopts::value<string>())
//        ("v,verbose", "Verbose output.", cxxopts::value<bool>()->default_value("false"))
        ("h,help", "Print usage")
    ;

    auto opts_parsed = opts.parse(argc, argv);

    if (argc == 1 || opts_parsed.count("help")){
        std::cerr << opts.help() << std::endl;
        return 1;
    }


    bool to_stdout = false;
    string seq_file = opts_parsed["i"].as<string>();
    string barcode_file = opts_parsed["b"].as<string>();
    string output_file;
    try{
        output_file = opts_parsed["o"].as<string>();
    } catch(cxxopts::exceptions::option_has_no_value& e){
        to_stdout = true;
    }

    if(to_stdout) filter_barcodes(seq_file, barcode_file, cout);
    else{
        ofstream out(output_file);
        filter_barcodes(seq_file, barcode_file, out);
    }

    return 0;
}


int main(int argc, char** argv){

    vector<string> commands = {"analyze", "filter"};
    if(argc == 1 || argv[1] == string("--help") || argv[1] == string("-h")){
        cerr << "Available commands: " << endl;
        for(string S : commands) cerr << "   " << argv[0] << " " << S << endl;
        cerr << "Running a command without arguments prints the usage instructions for the command." << endl;
        return 1;
    }

    string command = argv[1];

    // Drop the first element of argv
    for(int64_t i = 1; i < argc; i++) argv[i-1] = argv[i];
    argc--;

    if(command == "analyze") analyze_main(argc, argv);
    else if(command == "filter") filter_main(argc, argv);
    else{
        cerr << "Invalid command " << command << endl;
        return 1;
    }

    return 0;
}