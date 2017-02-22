//
//  csvoutput.cpp
//  searchNoise
//
//  Created by carushi on 11/09/22.
//  Copyright 2011年 Undergraduate Student, Department of Science, University of Tokyo. All rights reserved.
//

#ifndef searchNoise_csvoutput_cpp
#define searchNoise_csvoutput_cpp

#include "filename.h"


namespace mousebrain
{
    void Paint::csvoutput(const char* filename, class SmoothExpression &filedata)
    {
        string outputfile = RawFilename::numToCSV(filename);
        ofstream ofs(outputfile.c_str());
        ofs << "x, y, z, r, g, b" << endl;

        for (int i = 0; i < Filename::XNUM; i++)
            for (int j = 0; j < Filename::YNUM; j++)
                for (int k = 0; k < Filename::ZNUM; k++) {
                    double r,g,b;
                    if ( filedata.getData( i, j, k ) != 0.0 )
                    {
//                    setRGB( &r, &g, &b, filedata.getData(i, j, k), filedata.getMax());
                        ofs << i << ',' << j << ',' << k << ',' << r << ',' << g << ',' << b << std::endl;
                    }
                }
        return;
    }
}

#endif
