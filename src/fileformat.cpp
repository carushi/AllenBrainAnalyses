//
//  fileformat.cpp
//  searchNoise
//
//  Created by carushi on 11/10/27.
//  Copyright 2011年 Undergraduate Student, Department of Science, University of Tokyo. All rights reserved.
//

#include "fileformat.h"

namespace mousebrain
{
    SmoothExpression::~SmoothExpression() { freeData(); }
    void SmoothExpression::freeData()
    {
        if ( data == 0 ) return;
        for (int i = 0; i < XNUM; i++) {
            for (int j = 0; j < YNUM; j++)
                delete[] data[i][j];
            delete[] data[i];
        }
        delete[] data;
        return;
    }

    void SmoothExpression::clearData(double init)
    {
        if ( data == 0 ) return;
        for (int i = 0; i < XNUM; i++)
            for (int j = 0; j < YNUM; j++)
                for (int k = 0; k < ZNUM; k++)
                    data[i][j][k] = init;
        return;
    }

    void SmoothExpression::reserveData(double init)
    {
        data = new double**[XNUM];
        for (int i = 0; i < XNUM; i++) {
            data[i] = new double*[YNUM];
            for (int j = 0; j < YNUM; j++)
            {
                data[i][j] = new double[ZNUM];
                for (int k = 0; k < ZNUM; k++)
                    data[i][j][k] = init;
            }

        }
        return;
    }

    void SmoothExpression::printData()
    {
        for (int i = 0; i < XNUM; i++) {
            for (int j = 0; j < YNUM; j++) {
                for (int k = 0; k < ZNUM; k++)
                    cout << getData(i, j, k) << " ";
                cout << endl;
            }
            cout << endl;
        }
        return;
    }

    int SmoothExpression::getFileData(const char* filename)
    {
        int previewz = 0, count = 0, maxcount = 0;
        string str = numToStr(filename);
        ifstream ifs(str.c_str());
        getline(ifs,str);      getline(ifs,str);    /* header */
        if ( str != RawFilename::dimentionStr(XNUM, YNUM, ZNUM) ) return 1;
        while( getline(ifs, str) ) {
            Words v = myParse( str, "," );
            double c = setData(atoi(v[0].c_str()), atoi(v[1].c_str()), atoi(v[2].c_str()), atof(v[3].c_str()));
            if ( atoi(v[2].c_str()) != previewz ) {
                previewz = atoi(v[2].c_str());
                count = 0;
            }
            if ( ++count > maxcount ) {
                maxZ = previewz;
                maxcount = count;
            }
            if ( c > max ) max = c;
            if ( c > 0.0 && ( c < min || min == 0.0 )) min = c;
        }
        return 0;
    }

    SmoothExpressionZ::SmoothExpressionZ(int X, int Y, int Z) : SmoothExpression(X, Y, Z) {}

    void SmoothExpressionZ::setParam(int num, int &countmax, int &imax, int &jmax)
    {
        if ( num == ZNUM ) { countmax = XNUM*YNUM/100; imax = XNUM; jmax = YNUM; }
        else if ( num == XNUM ) { countmax = YNUM*ZNUM/100; imax = YNUM; jmax = ZNUM; }
        else { countmax = ZNUM*XNUM/100; imax = ZNUM; jmax = XNUM; }
        return;
    }

    double SmoothExpressionZ::getDataOpt(int num, int i, int j, int k)
    {
        if ( num == ZNUM ) return getData(i, j, k);
        else if ( num == XNUM ) return getData(k, i, j);
        else return getData(j, k, i);
    }

    void SmoothExpressionZ::findStart(int num) /* num means select type(x,y,z) and countmax means min space*/
    {
        int tmps = -1, tmpe = -1, count = 0, countmax, imax, jmax;
        setParam(num, countmax, imax, jmax);
        for ( int k = 0; k < num; k++, count = 0 ) {
            for ( int j = 0; j < jmax && count < countmax; j++ )
                for ( int i = 0; i < imax && count < countmax; i++ )
                    if ( getDataOpt(num, i, j, k) > 0 ) count++;
            if ( count < countmax ) continue;
            if ( tmps < 0 ) tmps = k;
            tmpe = k;
        }
        start = tmps;
        end = tmpe;
        return;
    }
}

