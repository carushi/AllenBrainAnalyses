//
//  filename.cpp
//  searchNoise
//
//  Created by carushi on 11/10/27.
//  Copyright 2011年 Undergraduate Student, Department of Science, University of Tokyo. All rights reserved.
//


#include "filename.h"

namespace mousebrain
{
    const double Filename::GAMMA = 0.2;
    const double Filename::PVALUE = 0.05;
    const bool Filename::init = false;
    const bool Filename::writeseedeffect = true;
    const bool Filename::miRNA = false;

    void Filename::readFilename(int option, int second, int third)
    {
        open(RawFilename::filename.c_str());
        switch(option) {
            case XOPTION:
                searchDataRange(XNUM); break;
            case YOPTION:
                searchDataRange(YNUM); break;
            case ZOPTION:
                searchDataRange(ZNUM); break;
            case SAMEAREA:
            case ALLAREA:
            case ALLFORMIRNA:
                readFilenameSameArea(option, second, third); break;
            case NORMAL: readFilenameNormal();
        }
        return;
    }

    void Filename::readFilenameSameArea(int option, int second, int third)
    {
        IDZLIST maxZ;
        class SmoothExpression tmpdata(XNUM, YNUM, ZNUM);
        getFilenameFromText(tmpdata, maxZ);
        if ( maxZ.size() == 0 ) return;
        cout << "(datasize = " << maxZ.size() << ")" << endl;

        reopen(RawFilename::filename.c_str());
        switch (option) {
                //case SAMEAREA:
                //    chooseArea(maxZ); break;
            case ALLAREA:
                analyzeCorrelation(maxZ, second, third); break;
            case ALLFORMIRNA:
                IDLIST filenameList;
                getFilenameList(maxZ, filenameList);
                analyzeMiRNAReg(filenameList, second, third);
        }
        return;
    }

    void Filename::searchDataRange(int maxnum) /* get start point, end point and length for each gene */
    {
        int *start = new int[maxnum], *end = new int[maxnum], *length = new int[maxnum];
        memset(start, 0, sizeof(int)*maxnum); memset(end, 0, sizeof(int)*maxnum);
        memset(length, 0, sizeof(int)*maxnum);
        class SmoothExpressionZ tmpdata(XNUM, YNUM, ZNUM);
        Words v;
        while ( fileRead(ifs, v, " \t") )
            for ( Words::iterator it = v.begin(); it != v.end(); it++ )
            {
                if ( tmpdata.getFileData((*it).c_str()) == 0 ) {
                    tmpdata.findStart(maxnum);
                    int tmps = tmpdata.getStart(), tmpe = tmpdata.getEnd();
                    if ( tmps >= 0 && tmps < maxnum ) start[tmps]++;
                    if ( tmpe >= 0 && tmpe < maxnum ) {
                        length[tmpe-tmps]++;
                        end[tmpe]++;
                    }
                    tmpdata.clearData(-1.0);
                }

            }
        cout << "s e l" << endl;
        for (int i = 0; i < maxnum; i++)
            cout << start[i] << " " << end[i] << " " << length[i] << endl;
        delete[] start; delete [] end; delete[] length;
        return;
    }

    void Filename::analyzeCorrelation(IDZLIST &maxZ, int foldchange, int z)
    {
        class PlaneExpressionForCorrelation X(maxZ, GENEMAX);
        X.reserveData();
        cout << "-analyze correlation" << endl;
        if ( foldchange ) {
            //X.outputCorForClustering();
            Paint::psoutput("heatmapgamma=0.2", X, threshold);
            //X.outputCorWithRandomization(true, true);
            //X.writeFoldChangeData(GAMMA);
            //X.outputCorWithRandomization(false, true);

        } else {
            X.outputCorWithRandomization(true, false);
            X.outputCorWithRandomization(false, false);
        }
        X.freeData();
        return;
    }

    void Filename::analyzeMiRNAReg(IDLIST &filenameList, int getz, int miRNANo)
    {
        cout << "-analyzeMiRNAReg" << endl;
        class GeneSequence geneSequence(false, filenameList);
        if (init) setSource(filenameList, geneSequence);
        else if (miRNA) {
            //ofstream ofs("/Users/username/program/atlas/data/searchNoise/output/source/dataNum.txt");
            //ofs << "z\tx\ty\t\tnum" << endl;
            class PlaneExpressionForHGD X(filenameList, &geneSequence, GENEMAX);
            X.initRes(DATASTART, true);
            analyzePvalue(X);
            X.freeData();
        } else threadGetPvalue(filenameList, getz, miRNANo);

        return;
    }

    Matrix Filename::getPvalue(class PlaneExpressionForHGD& X, int seed, int shfseq, bool up)
    {

        Matrix mindata = vector<vector<double> >(Filename::XNUM, vector<double> (Filename::YNUM, 1.0));
        for (int x = 0; x < XNUM; x++)
            for (int y = 0; y < YNUM; y++)
                if ( X.isData(x, y) )
                    mindata[x][y] = X.getMinPvalue(x, y, up);
        return mindata;
    }
    void Filename::outputPvalue(int seed, int jobnum, bool miRNA, class PlaneExpressionForHGD& X, Matrix& dlow, Matrix& dup, Quantmtx* qlow, Quantmtx* qup, int shuffle)
    {
        string filename = RawFilename::getpvalueFile(X.getz, jobnum);
        X.mutexLock();
        ofstream ofs(filename.c_str(), ios::app);
        for (int i = 0; i < Filename::XNUM; i++)
            for (int j = 0; j < Filename::YNUM; j++)
                if ( X.isData(i, j) ) {
                    /*
                    ofs << seed << "\t" << miRNA << "\t" << "\t" << i << "\t" << j << "\t" << X.getz << "\t"
                    << dlow[i][j] << "\t" << dup[i][j] << "\t";
                    if (qlow != NULL || qup != NULL)
                        ofs << (*qlow)[i][j]/10.0 << "\t" << (*qup)[i][j]/10.0 << "\t" << shuffle << endl;
                    else
                        ofs << "0\t0\t" << shuffle << endl;
                     */
                    ofs << seed << "\t" << i << "\t" << j << "\t" << -log(dlow[i][j]) << "\t" << -log(dup[i][j]) << "\t" << shuffle << endl;
                }
        ofs.close();
        X.mutexUnlock();
    }

    void Filename::countQuantile(Quantmtx& quant, const Matrix& oridata, const Matrix shfdata)
    {
        for (int i = 0; i < Filename::XNUM; i++)
            for (int j = 0; j < Filename::YNUM; j++)
                if ( oridata[i][j] < shfdata[i][j] ) quant[i][j]++;
    }

    void Filename::plotPoint(class PlaneExpressionForHGD &X, int threadnum, int jobnum)
    {
        //vector<pair<int, int> > seedList = readRandomSeed();
        vector<int> seedList = readAllSeed();
        vector<int>::iterator it = seedList.begin();
        for (int i = 0; i < MAXSEED; i++) {
            if ( i%MAXJOB != jobnum || (i/MAXJOB)%THREADNUM != threadnum ) continue;
            int seed = i, miRNAor = 0;
            for (; it != seedList.end() && *it <= seed; it++)
                if (*it == seed) miRNAor++;
            if (it != seedList.end()) cout << " nextseed : " << *it << " ";
            cout << "seed : " << i << "/" << MAXSEED << " " << miRNAor << endl;
            X.geneSequence->setSeedEffect(X.seedEffect, seed, 0, (miRNAor > 0));
            Matrix datalow = getPvalue(X, seed, 0, false), dataup = getPvalue(X, seed, 0, true);
            Quantmtx quantlow = vector<vector<int> >(Filename::XNUM, vector<int> (Filename::YNUM, 0)), quantup = quantlow;
            for (int j = 1; j < GeneSequence::MAXSHFSEQ+1 && j < Filename::SHFMAX; j++)
                {
                    X.geneSequence->setSeedEffect(X.seedEffect, seed, j, (miRNAor > 0));
                    Matrix shflow = getPvalue(X, seed, j, false), shfup = getPvalue(X, seed, j, true);
                    outputPvalue(seed, jobnum, miRNAor, X, shflow, shfup, NULL, NULL, j);
                    //countQuantile(quantlow, datalow, shflow);
                    //countQuantile(quantup, dataup, shfup);
                }
                outputPvalue(seed, jobnum, miRNAor, X, datalow, dataup, &quantlow, &quantup, 0);
                //outputPvalue(seed, jobnum, miRNAor, X, datalow, dataup, NULL, NULL, 0);
        }
        return;
    }

    void Filename::setSeed(class PlaneExpressionForHGD &X, int count, int shfseq)
    {
        X.geneSequence->setSeedEffect(X.seedEffect, count, shfseq, false);
        X.geneSequence->miRNATarget = X.geneSequence->getRandomSeed(count, GeneSequence::SEEDLENGTH);
        return;
    }

    void Filename::randomExtractOfSeed(ofstream& ofs, int randommax, vector<int>& noMiRNA)
    {
        if ( (int)noMiRNA.size() < randommax )
            for (vector<int>::iterator it = noMiRNA.begin(); it != noMiRNA.end(); it++)
                ofs << *it << "\t0" << endl;
        else {
            vector<int> alreadyCount((int)noMiRNA.size(), 0);
            mt19937 seedOfR((unsigned int)(time(NULL)));
            uniform_real<> range(0, (int)noMiRNA.size()-1);
            variate_generator<mt19937, uniform_real<> > rand(seedOfR, range);
            for (int count = 0; count < randommax; ) {
                int tempseed = (int)rand();
                if ( alreadyCount[tempseed] == 0 ) {
                    ofs << noMiRNA[tempseed] << "\t0" << endl;
                    alreadyCount[tempseed]++;
                    count++;
                }
            }
        }
    }

    void Filename::writeRandomSeed(class PlaneExpressionForHGD &X, int randommax)
    {
        ifstream ifs(RawFilename::miRNAfile.c_str());
        Words v, seed;
        vector<int> noMiRNA;
        fileRead(ifs, v, "\t");
        while ( ifs && fileRead(ifs, v, "\t") )
            if ( v.size() == 3 ) seed.push_back(v[2]);
        ofstream ofs(RawFilename::seedfile.c_str());
        ofs << "seed\tmiRNA" << endl;
        for (int count = 0; count < MAXSEED; count++) {
            string str = X.geneSequence->getRandomSeed(count, GeneSequence::SEEDLENGTH);
            if ( find(seed.begin(), seed.end(), str) != seed.end() )
                ofs << count << "\t1" << endl;
            else noMiRNA.push_back(count);
        }
        randomExtractOfSeed(ofs, randommax, noMiRNA);
        return;
    }

    vector<pair<int, int> > Filename::readRandomSeed()
    {
        ifstream ifs(RawFilename::seedfile.c_str());
        Words v;
        vector<pair<int, int> > seedList;
        fileRead(ifs, v, "\t");
        while (fileRead(ifs, v, "\t"))
            if ( v.size() == 2 )
                seedList.push_back(pair<int, int>(atoi(v[0].c_str()), atoi(v[1].c_str())));
        return seedList;
    }

    vector<int> Filename::readAllSeed()
    {
        ifstream ifs(RawFilename::seedfile.c_str());
        vector<int> seedBit;
        Words v;
        fileRead(ifs, v, "\t");
        while (fileRead(ifs, v, "\t"))
            if ( v.size() == 2 && atoi(v[1].c_str()) == 1 )
                seedBit.push_back(atoi(v[0].c_str()));
        sort(seedBit.begin(), seedBit.end());
        return seedBit;
    }


    void Filename::uniqGene()
    {
        cout << "--remove duplication" << endl;
        using namespace ushuffle;
        if ( RawFilename::local ) {
            system("/Users/username/Mylibrary/c/ushuffle/main > /Users/username/program/atlas/data/searchNoise/output/source/shuffle.txt");
        } else {
            class Unique unique;
            unique.uniq(RawFilename::filename7.c_str(), RawFilename::seqback.c_str());
            unique.uniqFasta(RawFilename::rawsequence.c_str(), RawFilename::rawseqback.c_str());
            //system("gcc -O3 -c -o /home/username/ushuffle/ushuffle.o /home/username/ushuffle/ushuffle.c");
            //system("gcc -o /home/username/ushuffle/main -O3 /home/username/ushuffle/main.c /home/username/ushuffle/ushuffle.o");
            //system("/home/username/ushuffle/main > /grid/username/shuffle.txt");
        }
    }

    void Filename::setSource(IDLIST &filenameList, GeneSequence &geneSequence)
    {
        class PlaneExpressionForHGD X(filenameList, &geneSequence, GENEMAX);
        X.initRes(DATASTART, true);
        if ( writeseedeffect ) {
            cout << "--writeSE\n\n" << endl;
            class PlaneExpressionForHGD Y(filenameList, &geneSequence, GENEMAX);
            Y.initRes(DATASTART, true);
            //outputSESample(X,Y);
            outputseedSESample(X,Y);
            Y.freeData();
        } else {
        //cout << "--read sequence" << endl;
        //geneSequence.threadWriteSequence();
            cout << "--search SE" << endl;
            geneSequence.threadWriteSeedEffect(filenameList);
            uniqGene();
            writeRandomSeed(X, 2000);
        }
        X.freeData();
        cout << "init!" << endl;
    }

    void Filename::outputSESample(PlaneExpressionForHGD& X, PlaneExpressionForHGD& Y)
    {
        cout << "#\tseed\tcount0\tcount1\tsize" << endl;
        cout << "*\tseed\tshuffle\tdifshf1\tdifraw1\tdifsum" << endl;
        for (int i = 0; i < MAXSEED; i++) {
            X.geneSequence->setSeedEffect(X.seedEffect, i, 0, false);
            int count1 = 0, count0 = 0;
            for (int j = 1; j < GeneSequence::MAXSHFSEQ+1; j++) {
                Y.geneSequence->setSeedEffect(Y.seedEffect, i, j, false);
                int difshf1 = 0, difraw1 = 0;
                for (map<int,int>::iterator it = X.seedEffect.begin(); it != X.seedEffect.end(); it++) {
                    if ( Y.seedEffect[it->first] > 0 && it->second == 0 ) difshf1++;
                    else if ( Y.seedEffect[it->first] == 0 && it->second > 0 ) difraw1++;
                    if (j == 1) {
                        if (it->second > 0) count1++;
                        else count0++;
                    }
                }
                cout << "*\t" << i << "\t" << j << "\t" << difshf1 << "\t" << difraw1 << "\t" << difshf1+difraw1 << endl;
            }
            cout << "#\t" << i << "\t" << count0 << "\t" << count1 << "\t" << X.seedEffect.size() << endl;
        }
    }
    void Filename::outputseedSESample(PlaneExpressionForHGD& X, PlaneExpressionForHGD& Y)
    {
        cout << "#\tseed\tcount0\tcount1\tsize" << endl;
        cout << "*\tseed\tshuffle\tdifshf1\tdifraw1\tdifsum" << endl;
        vector<int> dif1 = vector<int>(MAXSEED,0);
        for (int i = 0; i < MAXSEED; i++) {
            X.geneSequence->setSeedEffect(X.seedEffect, i, 0, false);
            int count1 = 0, count0 = 0;
            for (int j = 1; j < GeneSequence::MAXSHFSEQ+1; j++) {
                Y.geneSequence->setSeedEffect(Y.seedEffect, i, j, false);
                int difshf1 = 0, difraw1 = 0;
                for (map<int,int>::iterator it = X.seedEffect.begin(); it != X.seedEffect.end(); it++) {
                    if ( Y.seedEffect[it->first] > 0 && it->second == 0 ) difshf1++;
                    else if ( Y.seedEffect[it->first] == 0 && it->second > 0 ) difraw1++;
                    if (j == 1) {
                        if (it->second > 0) count1++;
                        else count0++;
                    }
                }
                cout << "*\t" << i << "\t" << j << "\t" << difshf1 << "\t" << difraw1 << "\t" << difshf1+difraw1 << endl;
            }
            cout << "#\t" << i << "\t" << count0 << "\t" << count1 << "\t" << X.seedEffect.size() << endl;
        }
    }

    void Filename::analyzePvalue(PlaneExpressionForHGD& X)
    {
        cout << "--miRNA" << endl;
        class MiRNA miRNA(X);
        //miRNA.compressData(THREADNUM, DATASTART, DATAEND, SHFMAX, XNUM, YNUM);
        miRNA.readCompressData(SHFMAX);
    }

    void Filename::threadGetPvalue(IDLIST& filenameList, int getz, int job)
    {
        cout << "--get pvalue about z= " << getz << endl;
        if ( getz < DATASTART || getz >= DATAEND ) return;
        string filename = RawFilename::getpvalueFile(getz, job-1);
        ofstream ofs(filename.c_str());
        //ofs << "seed\tmiRNA\tx\ty\tz\tlow10~4\tup10~4\tquantlow\tquantup\tshuffle" << endl;
        ofs << "seed\tx\ty\tlow10~4\tup10~4\tshuffle" << endl;
        ofs.close();
        ofs.open(RawFilename::debugfile.c_str(), ios::trunc); ofs.close();
        vector<fileThreadArg> thread_arg = getThreadArg(filenameList, THREADNUM, getz, job-1);
        thread::Thread<fileThreadArg> thread;
        thread.variableThread(thread_arg, THREADNUM, threadFile);
        freeThreadArg(thread_arg);
        return;
    }


    void Filename::getFilenameFromRawdata(class SmoothExpression& tmpdata, IDZLIST& maxZ)
    {
        reopen(RawFilename::filename.c_str());
        Words v;
        while ( fileRead(ifs, v, "\t") )
            for ( Words::iterator it = v.begin(); it != v.end(); it++ )
            {
                if ( tmpdata.getFileData((*it).c_str()) == 0 )
                    maxZ.push_back(pair<int, int>(atoi((*it).c_str()), tmpdata.getMaxZ()));
                tmpdata.clearData(-1.0);
            }
        return;
    }

    void Filename::getFilenameFromText(class SmoothExpression& tmpdata, IDZLIST &maxZ)
    {
        ifstream ifs3(RawFilename::maxZfile.c_str());
        cout << "-read " << RawFilename::maxZfile << endl;
        Words v;
        while( fileRead(ifs3, v, " \n\\\t ") )
            maxZ.push_back(pair<int, int>(atoi(v[0].c_str()), atoi(v[1].c_str())));
        return;
    }

    void Filename::getFilenameList(IDZLIST& maxZ, IDLIST& filenameList)
    {
        for(IDZLIST::iterator it = maxZ.begin(); it != maxZ.end(); it++)
            filenameList.push_back((*it).first);
        return;
    }


    void Filename::readFilenameNormal()
    {
        Words v;
        class SmoothExpression tmpdata(XNUM, YNUM, ZNUM);
        while ( fileRead(ifs, v, " \n\\\t ") )
            for ( Words::iterator it = v.begin(); it != v.end(); it++ )
                if ( tmpdata.getFileData((*it).c_str()) == 0 )
                    Paint::psoutput( (*it).c_str(), tmpdata );
                    //Paint::csvoutput( (*it).c_str(), tmpdata );
        return;
    }
    vector<class fileThreadArg> Filename::getThreadArg(IDLIST &filenameList, int max, int z, int job)
    {
        cout << "-make data for threads" << endl;
        vector<class fileThreadArg> thread_arg;
        for (int i = 0; i < max; i++)
        {
            class GeneSequence* geneSequence = new GeneSequence(false, filenameList);
            class PlaneExpressionForHGD* X = new PlaneExpressionForHGD(filenameList, geneSequence, GENEMAX);
            X->initResSaveMem(z, true);
           class fileThreadArg filethreadarg(i, job, X, this);
            thread_arg.push_back(filethreadarg);
        }
        return thread_arg;
    }

    void Filename::freeThreadArg(vector<fileThreadArg>& thread_arg)
    {
        for (vector<fileThreadArg>::iterator it = thread_arg.begin(); it != thread_arg.end(); it++)
        {
            it->X->freeData();
            delete it->X->geneSequence;
            delete it->X;
        }
        return;
    }

    static void* threadFile(void *arg)
    {
        fileThreadArg *Arg = (fileThreadArg*)arg;
        Filename *filename = Arg->filename;
        filename->plotPoint(*(Arg->X), Arg->threadnum, Arg->job);
        return NULL;
    }

}
