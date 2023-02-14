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

#include "base.h"
#include "seq.h"
#include "barcode.h"

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

    int adapter = 0;
    int adapter_spacer = 0;
    int adapter_barcode_spacer = 0;

    int adapter_rev = 0;
    int spacer_adapter_rev = 0;
    int spacer_barcode_adapter_rev = 0;

    int records = -1;
    std::unordered_map<barcode::Barcode, int> barcode_counts = {};

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
        

        if (align_forward.adapter)
            adapter++;
        if (align_forward.spacer)
            adapter_spacer++;
        if (align_forward.barcode)
            adapter_barcode_spacer++;

        if (align_reverse.adapter)
            adapter_rev++;
        if (align_reverse.spacer)
            spacer_adapter_rev++;
        if (align_reverse.barcode)
            spacer_barcode_adapter_rev++;

/*
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
*/

    }

    std::cout << "reads: " << records << "\n"
        << "..  adapter  ................... = " << adapter << "\n" 
        << "..  adapter  .......  spacer  .. = " << adapter_spacer << "\n"
        << "..  adapter  barcode  spacer  .. = " << adapter_barcode_spacer << "\n"
        << ".................... -adapter .. = " << adapter_rev << "\n"
        << ".. -spacer   ....... -adapter .. = " << spacer_adapter_rev << "\n"
        << ".. -spacer  -barcode -adapter .. = " << spacer_barcode_adapter_rev << "\n";

    gzclose(input_file);
    gzclose(output_file);
}

