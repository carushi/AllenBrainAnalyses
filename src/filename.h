//
//  filename.h
//  searchNoise
//
//  Created by carushi on 11/09/20.
//  Copyright 2011年 Undergraduate Student, Department of Science, University of Tokyo. All rights reserved.
//

#ifndef searchNoise_filename_h
#define searchNoise_filename_h

#include <cmath>
#include <unistd.h>
#include <limits>
#include <cstdlib>
#include <boost/random.hpp>
//#include "planefileformat.h"
#include "fileformat.h"
//#include "genesequence.h"
#include "miRNA.h"
#include "thread.h"
#include "uniq.h"


namespace mousebrain
{
    using boost::mt19937;
    using boost::uniform_real;
    using boost::variate_generator;

    typedef vector<pair<int, int> > IDZLIST;
    typedef vector<int> IDLIST;
    typedef vector<vector<int> > Quantmtx;
    class fileThreadArg;
    class Filename
    {
    private:
        static const bool init, writeseedeffect, miRNA;
        double threshold;
        ifstream ifs;
    public:
        Filename(double Threshold) : threshold(Threshold) {}
        ~Filename() {}
        //bool operator=();
        static const int NORMAL = 0, SAMEAREA = 1, XOPTION = 2, YOPTION = 3, ZOPTION = 4,
            ALLAREA = 5, ALLFORMIRNA = 6,NOMIRNA = 0, XNUM = 67, YNUM = 41, ZNUM = 58,
            GENEMAX = 21230, THREADNUM = 7, MAXJOB = 8, SHFMAX = 100;
        static const double GAMMA, PVALUE;
        void readFilename(int, int, int);
        void readFilenameSameArea(int, int, int);
        void searchDataRange(int);
        void analyzeCorrelation(IDZLIST&, int, int);
        void analyzeMiRNAReg(IDLIST&, int, int);
        Matrix getPvalue(class PlaneExpressionForHGD&, int, int, bool);
        void countQuantile(Quantmtx&, const Matrix&, const Matrix);
        void outputPvalue(int, int, bool, class PlaneExpressionForHGD&, Matrix&, Matrix&, Quantmtx*, Quantmtx*, int);
        void plotPoint(class PlaneExpressionForHGD &, int, int);
        void setSeed(class PlaneExpressionForHGD &X, int count, int shfseq);
        void randomExtractOfSeed(ofstream&, int, vector<int>&);
        void writeRandomSeed(class PlaneExpressionForHGD &, int);
        vector<pair<int, int> > readRandomSeed();
        vector<int> readAllSeed();
        void uniqGene();
        void setSource(IDLIST &, GeneSequence &);
        void outputSESample(PlaneExpressionForHGD&, PlaneExpressionForHGD&);
        void outputseedSESample(PlaneExpressionForHGD&, PlaneExpressionForHGD&);
        void analyzePvalue(PlaneExpressionForHGD&);
        void threadGetPvalue(IDLIST&, int, int);
        void getFilenameFromRawdata(class SmoothExpression&, IDZLIST&);
        void getFilenameFromText(class SmoothExpression&, IDZLIST&);
        void getFilenameList(IDZLIST&, IDLIST&);
        void readFilenameNormal();
        vector<fileThreadArg> getThreadArg(IDLIST&, int, int, int);
        void freeThreadArg(vector<fileThreadArg>&);
        void open(const char* filename) { ifs.open(filename); }
        void reopen(const char* filename) { ifs.close(); ifs.open(filename, ios::in); }
    };

    class fileThreadArg
    {
    public:
        int threadnum, job;
        PlaneExpressionForHGD* X;
        Filename* filename;
        fileThreadArg(int Threadnum, int Job, PlaneExpressionForHGD* tempX, Filename* FileName) : threadnum(Threadnum), job(Job), X(tempX), filename(FileName){}
        ~fileThreadArg(){}
        void setMutex(class Mutex* tempmutex) { X->setMutex(tempmutex); }
    };

    class Paint
    {
    public:
        Paint(){}
        ~Paint(){}
        static void setRGB(double*, double*, double*, double, double);
        static void setRGBmono(double*, double*, double*, double, double);
        static void setRGBheatmap(double*, double*, double*, double, double);
        static void setPvalue( double *r, double *g, double *b, double data, double bon, int ann );
        static void drawPs(ofstream &, int, int, int, double, double, double);
        static void setColor(int type, double &r, double &g, double &b, double bon, double data, int ann);
        static double getValue(int type, int x, int y, class PlaneExpressionForHGD &filedata);
        static void drawColorChart(ofstream &ofs, double r, double g, double b, int num);
        static void psoutput(const char*, class SmoothExpression &);
        static void psoutput(const char*, class PlaneExpressionForCorrelation &, double);
        static void psoutput(string, class PlaneExpressionForHGD &, int, int);
        static void csvoutput(const char*, class SmoothExpression &);
    };

    static void *threadFile(void *);

}
#endif
