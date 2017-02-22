//
//  getstring.h
//  searchNoise
//
//  Created by carushi on 11/10/28.
//  Copyright 2011年 Undergraduate Student, Department of Science, University of Tokyo. All rights reserved.
//

#ifndef searchNoise_getstring_h
#define searchNoise_getstring_h

#include <vector>
#include <sstream>
#include <fstream>
#include <boost/tokenizer.hpp>

namespace mousebrain
{
    using std::vector;
    using std::string;
    using std::ifstream;
    using std::ostringstream;
    typedef vector<string> Words;

    bool fileRead(ifstream &, Words &, const char*);
    bool fileReadWithoutSkip(ifstream &, Words &, const char*);
    Words myParse(string str, const char* sep_str);
    Words parseWithoutSkip(string str, const char* sep_str);
    class RawFilename
    {
    private:
    public:
        RawFilename() {}
        virtual ~RawFilename() {}
        static const int DIV;
        static const bool local, mask, up, wobblepair;
        static const char *rootdir, *chrdir, *csvdir, *cordir, *psdir;
        static string filename1, filename2, filename3, filename32, filename4, filename5, filename7, seqback,
        filename8, shuffle, filename14, filename16, filename17, filename18, filename19, filename20, structureFile, clusterFile,
        svadir, datafile, pvaluedir,filename, maxZfile, miRNAfile, seedfile,_annotationDataFileName,_clusteringDataFileName,
        _outputFilename, threebase, debugfile, pvalueFile, pvaluefile, sortfile, sortwfile, debugfile2, organfile, rawsequence, rawseqback, allSE;

        static void setFilename();
        static string getFilename(string, const char*);
        static string getFilename(const char*, const char*);
        static string countToStr(vector<int>&);
        static string numToBF(int z);
        static string numToSEBF(int);
        static string dimentionStr(int x, int y, int z);
        static string numToCSV(const char* filename);
        static string zToPvalue(int);
        static string zToFilename(int, bool);
        static string zToFilenameFC(int z, bool rank, bool foldchange);
        static string zgammaToFilename(int, double, bool);
        static string psFilename(const char*, double);
        static string faFilename(const char *, const char *);
        static string medianAndVariance(int);
        static string heatmapAndGamma(int, double);
        static string pvalue(string&, int);
        static string getpvalueFile(int, int);
        static string getpvalueSortFile(int);
        static string getTFile(int);
        static string getSeedFile(int);
        static void initSeedFile();
    };

    template <class File>
    string numToStr(const File filename)
    {
        ostringstream stream;
        stream << RawFilename::svadir << filename << ".sva";
        return stream.str();
    }


}

#endif
