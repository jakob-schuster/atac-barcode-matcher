#include <string>
#include <vector>

#include "ssw/ssw_cpp.h"

struct TotalResult
{
    bool adapter;
    bool spacer;
    bool barcode;
    std::string identifier;
};

namespace base
{
    using Base = char;
}

namespace seq
{
    // a sequence is just a string
    // (because it would be annoying to re-implement alignment algorithms)
    using Seq = std::string;
}

namespace barcode
{
    // a barcode is a sequence
    using Barcode = seq::Seq; 
    using Barcodes = std::vector<Barcode>;
}

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

namespace barcode
{
    const Barcode no_barcode = "";
    const size_t length = 16;

    struct Result
    {
        bool    found;
        int     dist;
        Barcode barcode;
    };

    Barcodes
    load_barcodes(const std::string filename);

    Barcodes
    rev_comp(const Barcodes& barcodes);

    struct Result
    align(const seq::Seq& seq, const Barcodes& barcodes, const double edit_dist_proportion);
};