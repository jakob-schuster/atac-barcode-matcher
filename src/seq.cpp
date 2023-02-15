#include <math.h>
#include <iostream>

#include "ssw/ssw_cpp.h"
#include "edlib/include/edlib.h"

#include "types.h"
#include "base.h"
#include "seq.h"

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

        std::cout << "rev comped " << seq << " into " << new_seq << "\n";
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
        std::cout << "seq remaining is " << seq << "\n";

        // find the spacer
        auto spacer_result = align(seq_remaining, spacer, edit_dist_proportion);
        if (!spacer_result.found)
            return total;
        
        total.spacer = true;

        auto spacer_location = strand == '+' ?
            adapter_result.start + spacer_result.start :
            spacer_result.start;
        std::cout << "\tspacer at " << spacer_location << " to " << spacer_result.end << "\n";
        std::cout << "\tbarcode subseq is from " << spacer_result.end << " to " << seq_remaining.size() - spacer_result.end << "\n";

        // cut down the sequence again
        seq_remaining = strand == '+' ?
            seq_remaining.substr(0, spacer_result.start) :
            seq_remaining.substr(spacer_result.end, seq_remaining.size() - spacer_result.end);
        
        std::cout << "\tbarcode subseq is " << seq_remaining.size() << "\n";
        std::cout << "seq remaining is " << seq_remaining << "\n";

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
        determines whether a fuzzy match of key is present in the first search_size characters of seq,
        using a simple implementation of Ukkonen

        DEPRECATED: NO LONGER USED because we align using edlib instead
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
