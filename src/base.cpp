#include "base.h"

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
