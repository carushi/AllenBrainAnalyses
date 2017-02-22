//
//  getrawsequence.h
//  searchNoise
//
//  Created by carushi on 11/11/29.
//  Copyright 2011年 Undergraduate Student, Department of Science, University of Tokyo. All rights reserved.
//

#ifndef searchNoise_getrawsequence_h
#define searchNoise_getrawsequence_h

#include "getstring.h"

namespace mousebrain
{
    using std::vector;
    using std::pair;
    typedef unsigned long ulong;
    typedef vector<pair<ulong, ulong> > ExonList;

    class Sequence
    {
    private:
        const char* filename;
    public:
        Sequence(){}
        ~Sequence(){}
        ExonList getExonList(ulong, ulong, ExonList&);
        string getFileName(const char*, string&);
        string getRawSequence(const char*, string&, ulong, ulong, ExonList &);
        ulong isMax(ulong a, ulong b, ulong c)
        {
            ulong max = ( a > b ) ? a : b;
            if ( max > c ) return max;
            else return c;
        }
        ulong isMin(ulong a, ulong b, ulong c)
        {
            ulong min = ( a < b ) ? a : b;
            if ( min < c ) return min;
            else return c;
        }
    };
}



#endif
