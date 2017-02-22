//
//  miRNA.h
//  searchNoise
//
//  Created by carushi on 12/01/08.
//  Copyright 2012年 Undergraduate Student, Department of Science, University of Tokyo. All rights reserved.
//

#ifndef searchNoise_miRNA_h
#define searchNoise_miRNA_h


#include "planefileformat.h"
#include "genesequence.h"
#include "mymutex.h"
#include "annotationtree.h"
#include "pvalue.h"

namespace mousebrain
{
    using mutex::Mutex;

    typedef vector<vector<double> > Matrix;
    typedef vector<vector<int> > NMatrix;
    typedef pair<double, int> Element;
    extern bool sort_first(const Element& first, const Element& second);
    bool sort_pval(const PvalPointData& first, const PvalPointData& second);

    class MiRNA : public RawFilename
    {
    private:
        vector<vector<int> > annotationList;
        class PlaneExpressionForHGD& X;
        static const double NULLHYP;
    public:
        explicit MiRNA(class PlaneExpressionForHGD&);
        ~MiRNA(){}
        static const double PERCENT;
        class Mutex *mutex;
        int getRank(double, vector<double>&);
        NMatrix getThreadAnnotationList(int);
        void getCount(int&, int&, int&, int, vector<pair<double, int> >&);
        void outputPercent(int, map<int, int>&, Matrix&, vector<vector<PvalData> >&, NMatrix&);
        void setBackGround(bool, int, map<int, int>&, Matrix&, vector<vector<PvalData> >&);
        void getTreeData(int);
        void writePercentThreadMemory(int);
        void writePercent();
        template<class Argument>
        void thread(int threadMax, vector<Argument>& thread_arg, void*(*function)(void*));
        void parsePval(int, vector<string>&, vector<PvalPointData>&);
        void outputAnnRate(int, Node*, vector<int>&, vector<int>&);
        void calcHGD(int, vector<PvalPointData>&, Node*);
        void readCompressData(int);
        void compressData(int, int, int, int, int, int);
        void mutexLock() { if ( mutex != NULL ) mutex->lock(); }
        void mutexUnlock() { if ( mutex != NULL ) mutex->unlock(); }
        bool isequal(int num, string& str) {
            if ( num == atoi(str.c_str()) ) return true;
            else return false;
        }
        bool range(int ann, int min, int max)
        {
            if ( ann >= min && ann <= max ) return true;
            else return false;
        }
    };


    class ThreadArgMiRNA
    {
    public:
        int threadnum;
        MiRNA* miRNA;
        ThreadArgMiRNA(int Threadnum, MiRNA* mirna) : threadnum(Threadnum), miRNA(mirna){}
        ~ThreadArgMiRNA(){}
        void setMutex(class mutex::Mutex *tempMutex) { miRNA->mutex = tempMutex; }
    };
    static void *threadMiRNA(void *);
}
#endif
