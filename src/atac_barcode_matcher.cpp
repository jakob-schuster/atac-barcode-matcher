#include "atac_barcode_matcher.h"
#include "ssw/ssw_cpp.h"
#include "edlib/include/edlib.h"

#include <string>
#include <fstream>
#include <algorithm>
#include <iostream>
#include <unordered_map>
#include <zlib.h>
#include <math.h>

#include "string.h"

namespace base
{
    /*
        takes the complement of one base in a Seq
    */
    Base
    comp(const Base& b) {
        switch (b) {
            case 'A':
                return 'T';
            case 'T':
                return 'A';
            case 'G':
                return 'C';
            case 'C':
                return 'G';
            default:
                return b;
        }
    }
};


namespace seq
{
    /*
        reverses a sequence,
        not in place
    */
    Seq
    rev_comp(const Seq& seq)
    {
        Seq new_seq = "";
        auto base = seq.end();

        for (int i = seq.size() - 1; i >= 0; i--) {
            new_seq += base::comp(seq[i]);
        }
        return new_seq;
    }


    /*
        gets the appropriate edit distance, given an edit distance proportion (typically 0.1) and a key
    */
    int
    get_edit_dist(const double edit_dist_proportion, const seq::Seq& key)
    {
        return ceil(key.size() * edit_dist_proportion);
    }

    struct Result
    align(const Seq& seq, const Seq& key, const double edit_dist_proportion)
    {
        int edit_dist = get_edit_dist(edit_dist_proportion, key);
        std::cout << "edit_dist = " << edit_dist << ", ";
        auto config = edlibNewAlignConfig(edit_dist, EDLIB_MODE_HW, EDLIB_TASK_LOC, NULL, 0);
        auto result = edlibAlign(key.c_str(), key.length(), seq.c_str(), seq.length(), config);
        std::cout << "editDistance = " << result.editDistance << "\n";

        // check if the key was not found in the seq
        if (result.editDistance == -1)
            return {false, 0, 0, 0};
        
        // otherwise, return the finalResult;
        struct Result final_result = {true, result.editDistance, result.startLocations[0], result.endLocations[0]};
        edlibFreeAlignResult(result);
        return final_result;
    }

    struct TotalResult
    align_strand(const Seq& seq, const barcode::Barcodes& barcodes, const char strand, const double edit_dist_proportion)
    {
        const Seq current_adapter = strand == '+' ?
            adapter :
            adapter_rev_comp;
        const Seq current_spacer = strand == '+' ? 
            spacer : 
            spacer_rev_comp;

        struct TotalResult total = {false, false, false, barcode::no_barcode};

        // find the adapter
        auto adapter_result = align(seq, current_adapter, edit_dist_proportion);
        if (!adapter_result.found)
            return total;
        
        total.adapter = true;
    
        std::cout << "\tadapter at " << adapter_result.start << "\n";

        // cut down the sequence
        auto seq_remaining = strand == '+' ?
            seq.substr(adapter_result.end, seq.size() - adapter_result.end) :
            seq.substr(0, adapter_result.start);
        seq_remaining = seq;

        // find the spacer
        auto spacer_result = align(seq_remaining, spacer, edit_dist_proportion);
        if (!spacer_result.found)
            return total;
        
        total.spacer = true;

        std::cout << "\tspacer at " << adapter_result.start + spacer_result.start << "\n";

        // cut down the sequence again
        seq_remaining = strand == '+' ?
            seq_remaining.substr(0, spacer_result.start) :
            seq_remaining.substr(spacer_result.end, seq_remaining.size() - spacer_result.end);
        seq_remaining = seq;
        
        std::cout << "\tbarcode subseq is " << seq_remaining.size() << "\n";
        // check that there's the right amount of sequence left
        if (abs((int)(seq_remaining.size() - barcode::length)) > 2)
            return total;
        
        // find the barcode
        auto barcode_result = barcode::align(seq, barcodes, edit_dist_proportion);
        if (!barcode_result.found)
            return total;
        
        total.barcode = true;
        total.identifier = barcode_result.barcode;
        
        std::cout << "\tbarcode is " << barcode_result.barcode << "\n";

        return total;
    }
   
    
    struct AdapterSpacer
    align_adapter_spacer(const Seq& seq, const char strand, const double edit_dist_proportion) 
    {
        const Seq current_adapter = strand == '+' ?
            adapter :
            adapter_rev_comp;
        const Seq current_spacer = strand == '+' ?
            spacer :
            spacer_rev_comp;
        
        auto adapter_result = align(seq, current_adapter, edit_dist_proportion);

        if (!adapter_result.found) {
            return {false, false};
        }

        auto seq_remaining = strand == '+' ?
            seq.substr(adapter_result.end, seq.size() - adapter_result.end) :
            seq.substr(0, adapter_result.start);
        
        auto spacer_result = align(seq_remaining, current_spacer, edit_dist_proportion);
        if (!spacer_result.found)
            return {true, false};

        return {true, true};
    }

    /*
        determines whether a fuzzy match of key is present in the first searchSize characters of seq,
        using a simple implementation of Ukkonen
    */
    bool
    contains_match(const Seq &seq, const Seq &key, const int search_size, const double &edit_dist_proportion)
    {
        struct Result
        result;

        int edit_dist = get_edit_dist(edit_dist_proportion, key);

        // preprocessing
        int m = key.length() - 1;
        int n = std::min((int)seq.length(), search_size) - 1;

        int C[m];
        for (int i = 0; i < m; ++i) {
            C[i] = i;
        }
        int lact = edit_dist + 1; // last active

        // searching
        for (int pos = 0; pos < n + 1; ++pos) {
            int Cp = 0, Cn = 0;
            for (int i = 0; i < lact + 1; ++i) {
                if (key[i] == seq[pos]) {
                    Cn = Cp;
                } else {
                    if (Cp < Cn) {
                        Cn = Cp;
                    }
                    if (C[i] < Cn) {
                        Cn = C[i];
                    }
                    Cn++;
                }
                Cp = C[i];
                C[i] = Cn;
            }

            // updating lact
            while (C[lact] > edit_dist) {
                lact--;
            }
            if (lact == m) {
                
                return 1;
            } else {
                lact++;
            }
        }

        return 0;
    }
}

namespace barcode
{
    /*
        loads the barcodes from a newline-separated file
    */
    Barcodes
    load_barcodes(const std::string filename) {
        std::ifstream
        file (filename);

        Barcodes
        barcodes = {};

        std::string barcode;
        while (std::getline(file, barcode)) {
            barcodes.push_back(barcode);
        }

        std::sort(barcodes.begin(), barcodes.end());
        return barcodes;
    }

    /*
        takes the reverse complememnt of all the barcodes
    */
    Barcodes
    rev_comp(const Barcodes& barcodes) {
        Barcodes
        rev_comp_barcodes = {};

        for (const auto& barcode : barcodes) {
            rev_comp_barcodes.push_back(seq::rev_comp(barcode));
        }

        return rev_comp_barcodes;
    }

    struct Result
    align(const seq::Seq& seq, const Barcodes& barcodes, const double edit_dist_proportion)
    {
        struct Result
        result = {false, 0, no_barcode};

        for (const auto & barcode : barcodes) {
            auto barcode_result = seq::align(seq, barcode, edit_dist_proportion);
            // auto barcode_found = seq::contains_match(seq, barcode, seq.size(), 0.1);
            if (barcode_result.dist != -1) {
                if (!result.found) {
                    result.found = true;
                    result.dist = 0;
                    result.barcode = barcode;
                }
            }
        }

        return result;
    }
};

void
read_line(
    gzFile file, char line [], int length
) {
    if (!gzgets(file, line, length)) {
        throw std::invalid_argument("Could not read line!");
    }
    line[strlen(line)-1]='\0';
}

void
write_record(
    gzFile file, char identifier [], char sequence [], char quality [], std::string barcode
) {
    int length = 500000;

    gzwrite(file, identifier, sizeof(char) * strlen(identifier));
    std::string new_identifier = " barcode=" + barcode;
    gzwrite(file, new_identifier.c_str(), sizeof(char) * new_identifier.length());
    gzwrite(file, std::string("\n").c_str(), sizeof(char) * 1);
    gzwrite(file, std::string("+\n").c_str(), sizeof(char) * 2);
    gzwrite(file, sequence, sizeof(char) * strlen(sequence));
    gzwrite(file, std::string("\n").c_str(), sizeof(char) * 1);
    gzwrite(file, quality, sizeof(char) * strlen(quality));
    gzwrite(file, std::string("\n").c_str(), sizeof(char) * 1);
}

void
write_line(
    gzFile file, char line [], int length
) {
    if (!gzgets(file, line, length)) {
        throw std::invalid_argument("Could not read line!");
    }
    line[strlen(line)-1]='\0';
}

int
main(int argc, char *argv[]) {
    const double edit_dist_proportion = 0.1;

    gzFile input_file = gzopen(argv[1], "r");
    gzFile output_file = gzopen(argv[2], "w");

    std::string barcodes_filename(argv[3]);

    auto barcodes = barcode::load_barcodes(barcodes_filename);
    auto rev_comp_barcodes = barcode::rev_comp(barcodes);

    int length = 500000;
    char identifier [length];
    char sequence   [length];
    char plus       [length];
    char quality    [length];

    int records = -1;
    int adapters = 0;
    int adapters_rev = 0;
    int spacers = 0;
    int spacers_rev = 0;
    int barcodes_count = 0;
    std::unordered_map<barcode::Barcode, int> barcode_counts = {};

    int adapter_spacer = 0;
    int spacer_adapter = 0;

    while (true) {
        std::cout << "record " << records << "\n";
        records++;

        try {
            read_line(input_file, identifier, length);
            read_line(input_file, sequence, length);
            read_line(input_file, plus, length);
            read_line(input_file, quality, length);
        }
        catch (const std::invalid_argument &e) {
            // if there's nothing to read, finish
            break;
        }

        std::string
        seq (sequence);
        std::cout << "record " << records << "\n";
        std::cout << seq << "\n";

        auto align_forward = seq::align_strand(seq, barcodes, '+', edit_dist_proportion);
        auto align_reverse = seq::align_strand(seq, barcodes, '-', edit_dist_proportion);
        
        auto adapter_spacer_forward = align_forward.adapter && align_forward.spacer;
        auto adapter_spacer_reverse = align_reverse.adapter && align_reverse.spacer;

        if (adapter_spacer_forward)
            adapter_spacer++;
        if (adapter_spacer_reverse)
            spacer_adapter++;

        auto valid_forward = adapter_spacer_forward && align_forward.barcode;
        auto valid_reverse = adapter_spacer_reverse && align_reverse.barcode;

        auto align_result = align_forward;
        if (valid_forward && !valid_reverse) {
            barcodes_count++;
            align_result = align_forward;
        } else if (!valid_forward && valid_reverse) {
            barcodes_count++;
            align_result = align_reverse;
        }

        if (valid_forward && valid_reverse)
            barcodes_count++;

        if (valid_forward || valid_reverse)
            write_record(output_file, identifier, sequence, quality, align_result.identifier);
        
        // auto forward = seq::align_simple(seq, alignment_stuff, '+');
        // auto reverse = seq::align_simple(seq, alignment_stuff, '-');
    }

    std::cout << "reads: " << records << "\n"
        << "..  adapter  .......  spacer  .. = " << adapter_spacer << "\n"
        << "..  adapter  barcode  spacer  .. = " << barcodes_count << "\n"
        << ".. -spacer   ....... -adapter .. = " << spacer_adapter << "\n"
        << ".. -spacer  -barcode -adapter .. = " << barcodes_count << "\n";

    gzclose(input_file);
    gzclose(output_file);
}

