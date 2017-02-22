//
//  genesequence.h
//  searchNoise
//
//  Created by carushi on 11/10/27.
//  Copyright 2011年 Undergraduate Student, Department of Science, University of Tokyo. All rights reserved.
//

#ifndef searchNoise_genesequence_h
#define searchNoise_genesequence_h

#include <algorithm>
#include <map>
#include <pthread.h>
#include <bitset>
#include "suffixarray.h"
#include "getrawsequence.h"
#include "baseString.h"
#include "thread.h"
#define MAXSEED (16384)

namespace mousebrain
{
    using std::string;
    using std::map;
    using std::ofstream;
    using std::ifstream;
    using std::cout;
    using std::endl;
    using std::ios;
    using std::fill;
    using std::bitset;
    using mutex::Mutex;
    using std::min;
    using std::max;

    typedef vector<int> IDLIST;
    typedef map<string, int> GeneList;
    typedef map<int, int> SeedEffect;
    typedef pair<int, int> EffectElement;
    typedef vector<pair<double, int> > Data;

    class MyThreadArg;
    class MyThreadSeedArg;
    class DNASeq
    {
    private:
    public:
        DNASeq(){}
        ~DNASeq(){}
        static const int ELEMENTNUM, SEEDLENGTH, THREAD, MAXMIRNA, MAXSHFSEQ;
        string getComplementarySequence(string&);
        void printNarray(int, string&, int *);
        string getSeedSequence(ifstream &, int, string&);
        string getSequence(Words&);
        bool istarget(string str, const char* target)
        {
            if ( strcmp(str.c_str(), target) == 0 ) return true;
            else return false;
        }

    };

    class GeneSequence : public DNASeq
    {
    private:
        class Mutex *mutex;
        vector<int> geneAnnotate;
        static const bool writemiRNA, debug;
    public:
        static const int MINLENGTH = 40, STOCKMAX = 100;

        GeneList  imageAndGaId;
        string miRNAName, miRNATarget;
        GeneSequence(bool newfilename, IDLIST& filenameList)
        {
            if ( newfilename ) writeFilename(filenameList);
            readFilename();
        }
        ~GeneSequence() {}
        int getCount(string&);
        string getRandomSeed(int, int);
        void writeSequence(int fileNum, string &sequence, const char* filename);
        void writeFileNum(int fileNum, const char* filename);
        void writeFilename(vector<int>&);
        bool stockCheck(const int, int, IDLIST*);
        void writeBitSE(int, vector<bitset<MAXSEED> >&);
        void readSeqFromTFile(int, const char*, IDLIST*);
        void readSeq(int threadnum, IDLIST*);
        void readFilename();
        void writeMiRNA(int, int);
        void* singleWrite(MyThreadSeedArg*);
        bitset<MAXSEED> countSeedExist(string, int);
        int* countThreeBaseOfSeed(string &);
        void countThreeBaseNum(string, int);
        void throwThread(IDLIST&);
        void threadWriteSeedEffect(IDLIST&);
        void readThreeBaseNum(Data& isid, long long int*);
        int getisid(string &str) { return atoi(str.substr(1).c_str()); }
        void setMutex(class Mutex* tempmutex) { mutex = tempmutex; }
        void mutexLock() { if ( mutex != NULL ) mutex->lock(); }
        void mutexUnlock() { if ( mutex != NULL ) mutex->unlock(); }
        const char* getDebugFile() { return RawFilename::filename3.c_str(); }

        void writeFileNumBinary(int isid, int shfseq)
        {
            ofstream ofs(RawFilename::numToSEBF(shfseq).c_str(), ios::binary|ios::app);
            ofs.write((char*)&isid, sizeof(int));
            ofs.close();
        }

        int getMin(int a, int b, int length)
        {
            if (a < 0) a = length;
            if (b < 0) b = length;
            return min(a, b);
        }

        void addSequence(vector<string>& tList, vector<string>& wList, string seq);
        vector<int> getWobble(int, bool, bool);

        inline void setSeedEffect(map<int, int> &seedEffect, int count, int shfseq, bool wobble)
        {
            seedEffect.clear();
            if ( count < 0 || count >= MAXSEED ) return;
            ifstream ifs(RawFilename::numToSEBF(shfseq).c_str(), ios::binary);
            vector<int> cList = getWobble(count, wobble, shfseq==0);
            bitset<MAXSEED> bs;
            for (int tempid; ifs.read((char*)&tempid, sizeof(int));  ) {
                if (ifs.read((char*)&bs, sizeof(bitset<MAXSEED>)) ) {
                    int exist = 0;
                    for (vector<int>::iterator it = cList.begin(); it != cList.end() && exist == 0; it++)
                        if (bs[*it] == 1 ) exist = 1;
                    seedEffect.insert(pair<int, int>(tempid, exist));
                }
            }
            return;
        }

    };



    class MyThreadArg
    {
    public:
        int threadnum;
        GeneSequence* geneSequence;
        class Mutex *mutex;
        MyThreadArg(int Threadnum, GeneSequence* GeneSeq) : threadnum(Threadnum), geneSequence(GeneSeq){}
        ~MyThreadArg(){}
        void setMutex(class Mutex *tempmutex) { geneSequence->setMutex(tempmutex); }
    };

    class MyThreadSeedArg : public MyThreadArg
    {
    public:
        IDLIST* filenameList;
        bool simple;
        bool writeMiRNA;
        MyThreadSeedArg(int Threadnum, GeneSequence* GeneSeq) :
            MyThreadArg(Threadnum, GeneSeq){}
        ~MyThreadSeedArg() {}
        void setFileList(IDLIST* filenamelist) { filenameList = filenamelist; }

    };

    static void *threadTrans(void *);
    static void callSystem(const char*, const char*, const char*);
    static void callSystem(const char* demand, const char* sourcefile);

}

#endif
