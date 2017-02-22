//
//  main.cpp
//  searchNoise
//
//  Created by carushi on 11/09/17.
//  Copyright 2011年 Undergraduate Student, Department of Science, University of Tokyo. All rights reserved.
//

#include <iostream>
//#include <cstdlib>
#include "filename.h"

#define SAMEAREA (1)
#define XOPTION (2)
#define YOPTION (3)
#define ZOPTION (4)
#define ALLAREA (5)
#define ALLFORMIRNA (6)

using namespace std;
using namespace mousebrain;


int main(int argc, const char * argv[])
{
    if ( argc < 2 ) return 1;
    cout << argv[0] << " " << argv[1] << " " << argv[2] << endl;
    switch (atoi(argv[1]))
    {
        case 1: { class Filename svafilename(0); svafilename.readFilename(XOPTION, 1, 1); return 0; }
        case 2: { class Filename svafilename(0); svafilename.readFilename(YOPTION, 1, 1); return 0; }
        case 3: { class Filename svafilename(0); svafilename.readFilename(ZOPTION, 1, 1); return 0; }
        case 4: if ( argc == 4 ) {
                    class Filename svafilename1(0.73);
                    svafilename1.readFilename(ALLAREA, atoi(argv[2]), atoi(argv[3]));
                    return 0;
                }
        case 5: if ( argc == 4 ) {
                    class Filename svafilename(0);
                    svafilename.readFilename(ALLFORMIRNA, atoi(argv[2]), atoi(argv[3]));
                return 0;
                }
        case 6: { class Filename svafilename(0); svafilename.readFilename(SAMEAREA, -1, 0); return 0; }
        default: return 0;
    }
}
