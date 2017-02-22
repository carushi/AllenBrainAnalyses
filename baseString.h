//
//  baseString.h
//  searchNoise
//
//  Created by carushi on 12/04/25.
//  Copyright 2012年 Undergraduate Student, Department of Science, University of Tokyo. All rights reserved.
//

#ifndef searchNoise_baseString_h
#define searchNoise_baseString_h

#include <string>
namespace mousebrain {

    using std::string;
    static void changeDNA(string&);
    static char basePair(char);
    static char getChar(int);
    static int getNum(char);
    static int encodeChar(char);

    static void changeDNA(string& str)
    {
        for (int i = 0; i < (int)str.length(); i++) {
            switch(str[i]) {
                case 'a': case 'A': str[i] = 'A'; break;
                case 'c': case 'C': str[i] = 'C'; break;
                case 'g': case 'G': str[i] = 'G'; break;
                case 't': case 'T': case 'u': case 'U': str[i] = 'T'; break;
                default: str[i] = 'N';
            }
        }
        return;
    }

    static char basePair(char c)
    {
        switch(c) {
            case 'a': case 'A': return 'T';
            case 'c': case 'C': return 'G';
            case 'g': case 'G': return 'C';
            case 't': case 'T': case 'u': case 'U': return 'A';
            default: return 'N';
        }
    }

    static char getChar(int count)
    {
        switch(count) {
            case 0: return 'A';
            case 1: return 'C';
            case 2: return 'G';
            case 3: default: return 'T';
        }
    }

    static int getNum(char c) /* 0~3 */
    {
        switch(c) {
            case 'A': case 'a': return 0;
            case 'C': case 'c': return 1;
            case 'G': case 'g': return 2;
            case 'T': case 't': return 3;
            default: return 0;
        }
    }
    static int encodeChar(char c)
    {
        switch ( c ) {
            case '$': return 0;
            case 'A': case 'a': case '1': return 1;
            case 'C': case 'c': case '2': return 2;
            case 'G': case 'g': case '3': return 3;
            case 'T': case 't': case '4': return 4;
            default: return 1;
        }
    }

}

#endif

