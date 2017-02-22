//
//  fileformat.h
//  searchNoise
//
//  Created by carushi on 11/09/17.
//  Copyright 2011年 Undergraduate Student, Department of Science, University of Tokyo. All rights reserved.
//

#ifndef searchNoise_fileformat_h
#define searchNoise_fileformat_h

#include <iostream>
#include <cstdlib>
#include "getstring.h"

namespace mousebrain
{
    using std::cout;
    using std::endl;

    class SmoothExpression
    {
    private:
        int maxZ;
        double max, min, ***data;
        double setData(int tmpx, int tmpy, int tmpz, double value)
        {
            if ( tmpx < 0 || tmpx >= XNUM || tmpy < 0 || tmpy >= YNUM || tmpz < 0 || tmpz >= ZNUM ) return -1.0;
            else  { data[tmpx][tmpy][tmpz] = value; return value; }
        }
    public:
        SmoothExpression(int X, int Y, int Z) : XNUM(X), YNUM(Y), ZNUM(Z){ reserveData(-1.0); max = min = 0.0; maxZ = 0; };
        ~SmoothExpression();
        const int XNUM, YNUM, ZNUM;
        int getMaxZ() { return maxZ; }
        int getFileData(const char* filename);
        double getMax() { return max; }
        double getData(int X, int Y, int Z)
        {
            if ( data == 0 || X < 0 || Y < 0 || Z < 0 || X >= XNUM || Y >= YNUM || Z >= ZNUM )
                return -1.0;
            else return data[X][Y][Z];
        }
        void initData() { data = 0; }
        void freeData();
        void reserveData(double init);
        void clearData(double init);
        void printData();
    };


    class SmoothExpressionZ : public SmoothExpression
    {
    private:
        int start, end;
    public:
        SmoothExpressionZ(int, int, int);
        ~SmoothExpressionZ(){};
        void print() { cout << start << " " << end << endl; }
        int getStart() { return start; }
        int getEnd() { return end; }
        void findStart(int num);
        void setParam(int num, int &countmax, int &imax, int &jmax);
        double getDataOpt(int num, int i, int j, int k);
    };


}

#endif
