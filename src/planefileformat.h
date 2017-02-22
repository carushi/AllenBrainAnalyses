//
//  planefileformat.h
//  searchNoise
//
//  Created by carushi on 11/10/11.
//  Copyright 2011年 Undergraduate Student, Department of Science, University of Tokyo. All rights reserved.
//

#ifndef searchNoise_planefileformat_h
#define searchNoise_planefileformat_h

#include <cstring>
#include <cmath>
#include <utility>
#include <limits>
#include <algorithm>
#include <ctime>
#include <boost/random.hpp>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/median.hpp>
#include <boost/accumulators/statistics/variance.hpp>
#include "annotationtree.h"
#include "genesequence.h"
#include "mymutex.h"

#define START (0.2)
#define END  (0.3)
#define STEP (1.0)
#define UPPER (1)
#define LOWER (2)
#define THREEBASE (64)
#define OUT (0.2)
#define DATASTART (6)
#define DATAEND (34)
#define DIAMETER (10.0)
#define HEATMAX (10.0)

namespace mousebrain
{
    using std::ofstream;
    using std::ifstream;
    using std::numeric_limits;
    using mutex::Mutex;
    typedef pair<double, int> Element;
    typedef vector<pair<int, int> > IDZLIST;
    typedef vector<vector<pair<int, int> > > Table;
    typedef map<int, int> SeedEffect;
    typedef vector<vector<int> > List;


    bool sort_first(const Element& first, const Element& second);
    bool sort_distance(const pair<double, double>& first, const pair<double, double>& second);
    bool sort_second(const Element& first, const Element& second);

    class Point
    {
    public:
        int x,y,z;
        Point(int X, int Y, int Z): x(X), y(Y), z(Z){}
        Point(){}
    };

    class PlaneExpression
    {
    private:
        double ***planeData;
    public:
        PlaneExpression(int GeneMax) : allGene(GeneMax){ dataGene = GeneMax; }
        virtual ~PlaneExpression(){}
        int getz, allGene, annMax, dataGene;
        static const int XNUM = 67, YNUM = 41, ZNUM = 58, MINPOINT = 1000, MINGENE = 212, MINGENE2 = 113;
        static const double THRESHOLD;
        List annotationData;
        void reserveData();
        void clearData();
        void freeData();
        void initDataSetting(int, bool, IDLIST*);
        void writeBF(bool);
        void printData() const;
        void printDataNum() const;
        void setAnnotation(int);
        void getRankVector(Data&, Data&, bool);
        void getMedianVariance(double***);
        void getPlaneDataFromBF(int, bool, IDLIST*);
        void setFoldChangeAgain(double, double***);

        const int getAnnData(int& x1, int& y1) const
        {
            if ( isRange(x1, y1) ) return annotationData[x1][y1]; else return 0;
        }
        bool isData(int& x1, int& y1) const
        {
            if ( isRange(x1, y1) && annotationData[x1][y1] != 0 )
                if (x1%2 == 0 && y1%2 == 0)
                    return true;
            return false;
        }
        bool isRange(int& x1, int& y1) const
        {
            if ( x1 >= 0 && x1 < XNUM && y1 >= 0 && y1 < YNUM ) return true;
            else return false;
        }
        bool isSameAnn(int ann1, int ann2) const
        {
            if ( ann1 == ann2 ) return true;
            else return false;
        }
        void setAllPlaneData(double ***newData) { planeData = newData; return; }
        double ***getAllPlaneData() { return planeData; }
        void setPointData(int gene, int X, int Y, double value)
        {
            if ( gene >= 0 && gene < dataGene && X >= 0 && X < XNUM && Y >= 0 && Y < YNUM )
                planeData[gene][X][Y] = value;
            return;
        }
        double distance(int x1, int y1, int x2, int y2)
        {
            return sqrt((double)((x1-x2)*(x1-x2))+(double)((y1-y2)*(y1-y2)));
        }
        double distance(int x1, int y1, int z1, int x2, int y2, int z2)
        {
            return sqrt((double)((x1-x2)*(x1-x2))+(double)((y1-y2)*(y1-y2))+(double)((z1-z2)*(z1-z2)));
        }
        double getPointData(int gene, int X, int Y) const
        {
            if ( gene >= 0 && gene < dataGene && X >= 0 && X < XNUM && Y >= 0 && Y < YNUM )
                return planeData[gene][X][Y];
            else return 0;
        }
        const int isSameAnn(int& x1, int& y1, int& x2, int& y2)
        {
            if (annotationData[x1][y1] != annotationData[x2][y2]) return 0;
            else return 1;
        }
        int getSameRank(int gene, Data& exprAndNum )
        {
            double value = exprAndNum[gene].first;
            for (int i = gene+1; i < (int)exprAndNum.size(); i++)
                if ( exprAndNum[i].first > value ) return i;
            return (int)exprAndNum.size();
        }
    };



    class PlaneExpressionForCorrelation : public PlaneExpression
    {
    private:
        IDZLIST& maxZlist;
        class Mutex *mutex;
    public:
        PlaneExpressionForCorrelation(IDZLIST& MaxZlist, int GeneMax) : maxZlist(MaxZlist), PlaneExpression(GeneMax) {}
        ~PlaneExpressionForCorrelation() {}
        void init(int, bool);
        boost::variate_generator<boost::mt19937, boost::uniform_real<> > getRand(int, int);
        void getOneGenePlaneData(int, int, int);
        void getPlaneData(int);
        void getDataForMaxZ();
        void setPointList(int, bool, vector<vector<pair<int, int> > > &);
        string printCorData(int, int, int, int, int, int, int, double);
        void outputAllCorrelation(bool, bool, int);
        void randomCalcCor(ofstream&, vector<Point>&, vector<Point>&, bool, bool, int);
        void allCalcCor(ofstream&, Table&, bool, bool);
        void annotationCalcCor(ofstream &, Table &, int, int, bool, bool);
        void outputCorWithRandomization(bool, bool);
        void outputCorForClustering();
        void outputCorWithAnnotate(bool, int, double);
        void appCorrelation(int, int, double, double, vector<pair<double, double> > &);
        void heatAndCorrelate(double);
        double calcMean(double, vector<pair<double, double> > &);
        double calcCor(Data&, Data&);
        double calcCorrelation(int, int, int, int, bool);
        double calcCorrelation(Data&, int, int, bool);
        double getHeatMapData(int, int, double);
        double getHeatMapDataDistance(int, int, double);
        void getMedianList(vector<double> &);
        void getMedianListNew(vector<double> &);
        void writeFoldChangeData(double);
        const char* getMVFile() { return RawFilename::filename18.c_str(); }
        const char* getCorFile(bool foldchange, bool rank) {
            if (foldchange)  {
                if (rank) return RawFilename::filename16.c_str();
                else return RawFilename::filename17.c_str();
            } else {
                if (rank) return RawFilename::filename19.c_str();
                else return RawFilename::filename20.c_str();
            }
        }
    };

    class PlaneExpressionForHGD : public PlaneExpression
    {
    private:
        class Mutex *mutex;
        double ***planeData;
    public:
        static const int BIN;
        bool reduce;
        GeneSequence *geneSequence;
        IDLIST &filenameList;
        map<int, int> seedEffect;
        PlaneExpressionForHGD(IDLIST &FilenameList, GeneSequence *Genesequence, int GeneMax) :
        filenameList(FilenameList), geneSequence(Genesequence), PlaneExpression(GeneMax) {
            mutex = NULL;
            reduce = false;
        }
        ~PlaneExpressionForHGD(){}
        int countData();
        void initRes(int, bool);
        void initResSaveMem(int, bool);
        void readFilenameHavSeq();
        void printDataNumOption() const;
        int threeBaseNum(string);
        int maxSeedGene(Data&, bool);
        double hyperGeographicDistribution(Data &, int, bool, bool);
        double getCoefficientOfM(SeedEffect &, string &, bool);
        void getCount(double *, double *, SeedEffect &);
        void getM(int &, int &, int, Data &, bool, bool);
        double getMinPvalue(int, int, bool);
        vector<double> outputPvalue(int, int);
        double lcombination(int N, int n) { return lgamma(N+1)-lgamma(n+1)-lgamma(N-n+1); }
        double probability(int N, int n, int m, int k)
        { return exp(lcombination(m,k)+lcombination(N-m,n-k)-lcombination(N,n)); }
        bool haveSequenceData(int gene)
        {
            if ( filenameList[gene] < 0 ) return false;
            else return true;
        }
        void mutexLock() { if ( mutex != NULL ) mutex->lock(); }
        void mutexUnlock() { if ( mutex != NULL ) mutex->unlock(); }
        void setMutex(class Mutex *tempMutex) { mutex = tempMutex; }
        const char* getDBFile() { return RawFilename::debugfile.c_str(); }

        inline double sumProbability(int N, int n, int m, int k, bool output)
        {
            double sum = 0.0;
            for ( ; k <= m && k <= n; k++) {
                double firstProbability = probability(N, n, m, k);
                if ( firstProbability == 0.0 ) continue;
                sum += firstProbability;
                for (k = k+1; k <= m && k <= n; k++) {
                    firstProbability = firstProbability*(double)(m-k+1)/(double)(k)*(double)(n-k+1)/(double)(N-m-n+k);
                    if ( firstProbability == 0.0 ) break;
                    sum += firstProbability;
                }
            }
            if ( output ) {
                mutexLock();
                ofstream ofs(RawFilename::debugfile.c_str(), ios::app);
                ofs << "N n m k" << " " <<  N << " " << n << " " << m << " "  << k << " " << sum << endl;
                mutexUnlock();
            }
            return sum;
        }


    };

}



#endif
