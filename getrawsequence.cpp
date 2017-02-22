//
//  getrawsequence.cpp
//  searchNoise
//
//  Created by carushi on 11/11/29.
//  Copyright 2011年 Undergraduate Student, Department of Science, University of Tokyo. All rights reserved.
//

#include "getrawsequence.h"
namespace mousebrain
{
    string Sequence::getFileName(const char* filename, string &chrom)
    {
        ostringstream stream;
        if (RawFilename::local)
            stream << filename << chrom << ".fa";
        else
            stream << filename << chrom << ".fa";

        return stream.str();
    }

    ExonList Sequence::getExonList(ulong start, ulong end, ExonList& exonSE)
    {
        ExonList exon;
        for (ExonList::iterator it = exonSE.begin(); it != exonSE.end(); it++)
            if ( it->first < it->second ) {
                if ( ( start <= it->first && it->first <= end ) || ( it->first < start && it->second >= start ) )
                    exon.push_back(std::pair<ulong, ulong>(it->first, it->second));
            } else {
                exon.clear();
                return exon;
            }
        return exon;
    }

    string Sequence::getRawSequence(const char* filename, string &chrom, ulong start, ulong end, ExonList& exonSE)
    {
        ulong count = 1;
        string str = getFileName(filename, chrom), sequence = "";
        ifstream ifs(str.c_str());
        ExonList exon = getExonList(start, end, exonSE);
        if ( !getline(ifs, str) ) return sequence;
        for (ExonList::iterator it = exon.begin(); it != exon.end() && getline(ifs, str); count += str.length()) {
            if ( count+str.length() < start || count+str.length() < it->first ) continue;
            ulong tmpstart = isMax(start, it->first, count),
            tmpend = isMin(end, it->second, count+str.length());
            sequence += str.substr(tmpstart-count, tmpend-tmpstart+1);
            if ( count+str.length() > it->second || count+str.length() > end ) it++;
        }
        return sequence;
    }
}
