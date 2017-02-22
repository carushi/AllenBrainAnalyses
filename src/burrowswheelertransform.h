//
//  burrowswheelertransform.h
//  searchNoise
///Users/username/Downloads/workspace/BWT/src/BWT.cpp
//  Created by carushi on 11/12/02.
//  Copyright 2011年 Undergraduate Student, Department of Science, University of Tokyo. All rights reserved.
//

#ifndef searchNoise_burrowswheelertransform_h
#define searchNoise_burrowswheelertransform_h

//#include <vector>
//#include <string>
#include <iostream>
#include "getstring.h"

#define OCCSIZE (10)
#define THREEBASE (64)

namespace mousebrain
{
    typedef vector<vector<int> > Array;
    using std::cout;
    using std::endl;
    using std::fill;

    class BWT {
    private:
        const int *s;
        const int *SA;
        const int slen;
        const int K;
        int *bkt;
    public:
        BWT(int *S, int *tempSA, int Slen, int tempK) : s(S), SA(tempSA), slen(Slen), K(tempK-1)
        { bkt = new int[K]; }
        ~BWT() { delete [] bkt; }
        int getOcc(int, int, Array&);
        int search(int*, int, Array&, bool);
        int BWTransform(string&);
        void getSeedNum(vector<string>&, vector<int>&);
        void printDEBUG(Array &, int, int, int*, int);
        void getQuery(string &str, int*);
        void getQueryThreeBase(int*, int);
        void makeOcc(Array&);
        int* getThreeBaseCount();
    };
}



#endif
