//
//  getstring.cpp
//  searchNoise
//
//  Created by carushi on 11/10/15.
//  Copyright 2011年 Undergraduate Student, Department of Science, University of Tokyo. All rights reserved.
//


#include "getstring.h"
namespace mousebrain
{
    const int RawFilename::DIV = 5;
    const bool RawFilename::local = false;
    const bool RawFilename::mask = false;
    const bool RawFilename::up = true;
    const bool RawFilename::wobblepair = false;

    const char* RawFilename::rootdir = (RawFilename::local) ?  "/Users/username/program/atlas/data/searchNoise/"
    :  "/home/username/rawdata/searchNoise/";
    string RawFilename::filename1 = getFilename(rootdir,"source/gene-series-map.tsv");
    string RawFilename::filename2 = getFilename(rootdir, "source/Refsec-gene.txt");
    string RawFilename::filename3 = getFilename(rootdir, "output/debug/shuffledebug.txt");
    string RawFilename::filename32 = getFilename(rootdir, "output/debug/normaldebug.txt");
    string RawFilename::filename4 = getFilename(rootdir, "source/miRNA.dat");
    string RawFilename::filename5 = getFilename(rootdir, "source/geneid.dat");
    string RawFilename::filename7 = getFilename(rootdir, "output/source/havesequence.txt");
    string RawFilename::seqback = getFilename(rootdir, "output/source/seqback.txt");
    string RawFilename::allSE = getFilename(rootdir, "output/source/allSE/se");
    string RawFilename::shuffle = "/grid/username/shuffle.txt";
    string RawFilename::filename16 = getFilename(rootdir, "output/correlation/rankcorrelation.txt");
    string RawFilename::filename17 = getFilename(rootdir, "output/correlation/expcorrelation.txt");
    string RawFilename::filename18 = getFilename(rootdir, "output/source/medianAndVariance.txt");
    string RawFilename::filename19 = getFilename(rootdir, "output/correlation/rankcorrelationrawdata.txt");
    string RawFilename::filename20 = getFilename(rootdir, "output/correlation/expcorrelationrawdata.txt");
    string RawFilename::threebase = (RawFilename::mask) ?
    getFilename(rootdir, "output/source/sevenbasethreebase.txt") : getFilename(rootdir, "output/source/threebase");
    string RawFilename::debugfile = (RawFilename::mask) ?
    getFilename(rootdir, "output/debug/pvaluesyldebug.txt") : getFilename(rootdir, "output/debug/pvaluealldebug.txt");
    string RawFilename::pvalueFile = (RawFilename::mask) ?
    getFilename(rootdir, "ourput/pvalue/pvaluesyllow.txt") : getFilename(rootdir, "output/pvalue/pvalueall.txt");
    string RawFilename::svadir = (RawFilename::local) ?
    "/Users/username/program/atlas/data/rawdata/" : "/home/username/rawdata/sva/";
    string RawFilename::datafile = getFilename(rootdir, "output/expression/");
    const char* RawFilename::csvdir = "output/csv/";
    const char* RawFilename::cordir = "output/correlation/";
    const char* RawFilename::psdir = "output/expression/";
    string RawFilename::pvaluedir = getFilename(rootdir, "output/pvalue/");
    string RawFilename::filename = getFilename(rootdir, "output/source/saggital_file.txt");
    string RawFilename::maxZfile = getFilename(rootdir, "output/source/maxZ2.txt");
    string RawFilename::miRNAfile = getFilename(rootdir, "output/source/miRNA.txt");
    string RawFilename::seedfile = getFilename(rootdir, "output/source/seed.txt");
    string RawFilename::structureFile = getFilename(rootdir, "source/brainstructure.csv");
    string RawFilename::clusterFile = getFilename(rootdir, "output/correlation/rankcorrelationforclustering.csv");
    const char* RawFilename::chrdir = (RawFilename::local) ?
    "/Users/username/program/data/atlas/chr/" : "/home/username/rawdata/chr/";
    string RawFilename::_annotationDataFileName = getFilename(rootdir, "source/AtlasAnnotation200.sva");
    string RawFilename::_clusteringDataFileName = getFilename(rootdir, "output/correlation/rankcorrelationwith2.0.dat");
    string RawFilename::_outputFilename = getFilename(rootdir, "output/debug/annotationoutput.txt");
    string RawFilename::pvaluefile = (up) ?
    ((wobblepair) ? getFilename(rootdir, "output/pvalue/wpvalueFile.txt") : getFilename(rootdir, "output/pvalue/pvalueFile.txt")) :
    ((wobblepair) ? getFilename(rootdir, "output/pvalue/wpvalueFilelow.txt") :getFilename(rootdir, "output/pvalue/pvalueFilelow.txt"));
    string RawFilename::sortfile = (up) ?
    getFilename(rootdir, "output/pvalue/pvalueFileSort.txt") : getFilename(rootdir, "output/pvalue/pvalueFilelowSort.txt");
    string RawFilename::sortwfile = (up) ?
    getFilename(rootdir, "output/pvalue/wpvalueFileSort.txt") : getFilename(rootdir, "output/pvalue/wpvalueFilelowSort.txt");
    string RawFilename::debugfile2 = (up) ?
    getFilename(rootdir, "output/debug/debug.txt") : getFilename(rootdir, "output/debug/debuglow.txt");
    string RawFilename::organfile = (wobblepair) ?
    getFilename(rootdir, "output/quantile/worgan.txt") : getFilename(rootdir, "output/quantile/organ.txt");
    string RawFilename::rawsequence = getFilename(rootdir, "output/source/rawsequence.txt");
    string RawFilename::rawseqback = getFilename(rootdir, "output/source/rawseqback.txt");


    string RawFilename::getFilename(string str, const char* filename)
    {
        string newfilename = str + filename;
        return newfilename;
    }
    string RawFilename::getFilename(const char* str, const char* filename)
    {
        string newfilename = str;
        newfilename += filename;
        return newfilename;
    }

    string RawFilename::countToStr(vector<int>& count)
    {
        ostringstream stream;
        for (vector<int>::iterator it = count.begin(); it != count.end(); it++)
            stream << *it << " ";
        return stream.str();
    }

    string RawFilename::numToBF(int z)
    {
        ostringstream stream;
        stream << datafile << "dataz" << z << "plane";
        return stream.str();
    }

    string RawFilename::numToSEBF(int shf)
    {
        ostringstream stream;
        stream << allSE << shf;
        return stream.str();
    }

    string RawFilename::dimentionStr(int x, int y, int z)
    {
        ostringstream stream;
        stream << "Dimensions:" << x << "," << y << "," << z;
        return stream.str();
    }

    string RawFilename::numToCSV(const char* filename)
    {
        ostringstream stream;
        stream << rootdir << csvdir << filename << ".csv";
        return stream.str();
    }

    string RawFilename::zToFilename(int z, bool rank)
    {
        ostringstream stream;
        stream << rootdir << cordir;
        if ( rank ) stream << "rankcorrelation";
        else stream << "expcorrelation";
        stream << z << ".txt";
        return stream.str();
    }

    string RawFilename::zToFilenameFC(int z, bool rank, bool foldchange)
    {
        ostringstream stream;
        stream << rootdir << cordir;
        if (rank) stream << "rankcorrelation";
        else stream << "expcorrelation";
        stream << z;
        if (!foldchange) stream << "withoutfoldchange";
        stream << ".txt";
        return stream.str();
    }

    string RawFilename::zToPvalue(int z)
    {
        ostringstream stream;
        stream << "Pvalue" << z;
        return stream.str();
    }

    string RawFilename::zgammaToFilename(int z, double gamma, bool rank)
    {
        ostringstream stream;
        stream << rootdir << cordir;
        if (rank) stream << "rankcorrelation";
        else stream << "expcorrelation";
        stream << z << "gamma" << gamma << ".txt";
        return stream.str();
    }

    string RawFilename::medianAndVariance(int z)
    {
        ostringstream stream;
        stream << rootdir << cordir << "medianAndVariance" << z << ".txt";
        return stream.str();
    }

    string RawFilename::psFilename(const char* filename, double correlation_min)
    {
        ostringstream stream;
        stream << rootdir << psdir << filename << correlation_min << ".ps";
        return stream.str();
    }

    string RawFilename::faFilename(const char *filename, const char *chromosome)
    {
        ostringstream stream;
        stream << filename << chromosome << ".fa";
        return stream.str();
    }

    string RawFilename::heatmapAndGamma(int getz, double gamma)
    {
        ostringstream stream;
        stream << "heatmapwithmonogamma" << gamma << "z" << getz;
        return stream.str();
    }

    string RawFilename::pvalue(string &miRNAName, int getz)
    {
        ostringstream stream;
        stream << "pvalue" << getz << miRNAName;
        return stream.str();
    }

    string RawFilename::getpvalueFile(int z, int jobnum)
    {
        ostringstream stream;
        stream << pvaluedir;
        if (wobblepair) stream << "w";
        stream << "array" << z << "job" << jobnum << ".txt";
        return stream.str();
    }

    string RawFilename::getpvalueSortFile(int z)
    {
        ostringstream stream;
        stream << pvaluedir;
        if (wobblepair) stream << "w";
        stream << "sort" << z << ".txt";
        return stream.str();
    }
    string RawFilename::getTFile(int num)
    {
        ostringstream stream;
        stream << rootdir << "output/temp" << num << ".txt";
        return stream.str();
    }

    string RawFilename::getSeedFile(int seed)
    {
        ostringstream stream;
        if (wobblepair)
            stream << rootdir << "output/source/wseed" << seed%DIV;
        else
            stream << rootdir << "output/source/seed" << seed%DIV;
        return stream.str();
    }
    void RawFilename::initSeedFile()
    {
        std::ofstream ofs;
        for (int i = 0; i < DIV; i++) {
            ostringstream stream;
            stream << rootdir << "output/source/seed" << i;
            string filename = stream.str();
            ofs.open(filename.c_str());
            ofs.close();
        }
    }

    Words parseWithoutSkip(string str, const char* sep_str)
    {
        typedef boost::tokenizer<boost::char_separator<char> > token;
        boost::char_separator<char> sep(sep_str, "", boost::keep_empty_tokens);
        token tok(str, sep);
        vector<string> v;
        for (token::iterator it = tok.begin(); it != tok.end(); it++) v.push_back(*it);
        return v;
    }

    Words myParse(string str, const char* sep_str)
    {
        typedef boost::tokenizer<boost::char_separator<char> > token;
        boost::char_separator<char> sep(sep_str);
        token tok(str, sep);
        Words v;
        for (token::iterator it = tok.begin(); it != tok.end(); it++) v.push_back(*it);
        return v;
    }

    bool fileRead(ifstream &ifs, Words &v, const char* sep)
    {
        string str;
        if ( !getline(ifs, str) ) return false;
        v = myParse(str, sep);
        return true;
    }

    bool fileReadWithoutSkip(ifstream &ifs, Words &v, const char* sep)
    {
        string str;
        if ( !getline(ifs, str) ) return false;
        v = parseWithoutSkip(str, sep);
        return true;
    }

}

