#ifndef BARCODE_H
#define BARCODE_H

#include <vector>

#include "types.h"
#include "seq.h"

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

#endif