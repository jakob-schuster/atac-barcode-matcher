#ifndef SEQ_H
#define SEQ_H

#include <string>

#include "types.h"
#include "barcode.h"

namespace seq
{
    // store the adapter and spacer sequences
    const Seq adapter = "AATGATACGGCGACCACCGAGATCTACAC";
    const Seq spacer  = "CGCGTCTG";
    const Seq adapter_rev_comp = "GTGTAGATCTCGGTGGTCGCCGTATCATT";
    const Seq spacer_rev_comp  = "CAGACGCG";

    struct AdapterSpacer
    {
        bool adapter;
        bool spacer;
    };

    struct Result
    {
        bool found;      // whether it's been found
        int  dist;       // edit distance
        int  start, end; // the start and end position of the result
    };

    Seq
    rev_comp(const Seq& seq);

    struct Result
    align(const Seq& seq, const Seq& key, const double edit_dist_proportion);

    struct TotalResult
    align_strand(const Seq& seq, const barcode::Barcodes& barcodes, const char strand, const double edit_dist_proportion);
}


#endif
