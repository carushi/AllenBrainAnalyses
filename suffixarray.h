//
//  suffixarray.h
//  searchNoise
//
//  Created by carushi on 11/11/29.
//  Copyright 2011年 Undergraduate Student, Department of Science, University of Tokyo. All rights reserved.
//

#ifndef searchNoise_suffixarray_h
#define searchNoise_suffixarray_h
//#include <iostream>
//#include <string>
//#include "getstring.h"
#include "burrowswheelertransform.h"
#define BASE (5)

namespace mousebrain
{

    class SuffixArray
    {
    private:
        int *str, *SA, length;
    public:
        SuffixArray(){}
        ~SuffixArray(){}
        int countSequence(string&, string&);
        string countSeed(string&, vector<string>&);
        void countSeed(string&, vector<string>&, vector<int>&);
        void printSA(int);
        void prints(int*, int);
        void printDebug(std::string &);
        void reserveArray(int, std::string&);
        void freeArray();
        void tempSetting(int*, int*, int);
        void inducedSorting(int*, int*, int*, int, int, bool);
        void getBucket(int*, int*, int, int, bool);
        bool moveLMS(int*, int*, int*, int, int);
        void SAIS(int*, int, int);
        int encodeChar( char c ) {
            switch ( c ) {
                case '$': return 0;
                case 'A': case 'a': case '1': return 1;
                case 'C': case 'c': case '2': return 2;
                case 'G': case 'g': case '3': return 3;
                case 'T': case 't': case '4': return 4;
                default: return 1;
            }
        }
        void tset(int *t, int i, int b, int Length) { if (i >= 0 && i < Length) t[i] = b; return; }
        int tget(int *t, int i, int Length) { if (i >= 0 && i < Length) return t[i]; else return -1; }
        int isLMS(int *t, int i, int Length) { return ( ( i > 0 && tget(t, i, Length) == 1 && tget(t, i-1, Length) == 0 ) ? 1 : 0); }
        int* countThreeBase(string &);
    };
}

#endif
