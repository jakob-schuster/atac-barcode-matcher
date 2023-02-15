#include <fstream>
#include <algorithm>

#include "types.h"
#include "barcode.h"

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
            if (barcode_result.dist != -1) {
                if (!result.found) {
                    result.found = true;
                    result.dist = barcode_result.dist;
                    result.barcode = barcode;
                } else if (barcode_result.dist < result.dist) {
                    result.dist = barcode_result.dist;
                    result.barcode = barcode;
                }
            }
        }

        return result;
    }
};