#ifndef TYPES_H
#define TYPES_H

#include <string>
#include <vector>

struct TotalResult
{
    bool adapter;
    bool spacer;
    bool barcode;
    std::string identifier;
};

namespace base
{
    // a base is one character
    // (because it would be expensive to re-implement it as an enum)
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


#endif