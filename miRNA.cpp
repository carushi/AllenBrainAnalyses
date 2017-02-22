//
//  miRNA.cpp
//  searchNoise
//
//  Created by carushi on 12/01/08.
//  Copyright 2012年 Undergraduate Student, Department of Science, University of Tokyo. All rights reserved.
//


#include "miRNA.h"

namespace mousebrain
{
    MiRNA::MiRNA(class PlaneExpressionForHGD& x): X(x) {}
    const double MiRNA::PERCENT = 0.95;
    const double MiRNA::NULLHYP = 0.05;

    bool sort_pval(const PvalPointData& first, const PvalPointData& second) {
        return (first.pvalue > second.pvalue);
    }

    int MiRNA::getRank(double pvalue, vector<double>& pvalueList)
    {
        for (vector<double>::iterator it = pvalueList.begin(); it != pvalueList.end(); it++)
            if ( pvalue > *it ) return (int)distance(pvalueList.begin(), it);
        return (int)pvalueList.size();
    }

    NMatrix MiRNA::getThreadAnnotationList(int num)
    {
        mutexLock();
        X.initDataSetting(num, true, &(X.filenameList));
        NMatrix annotationList = X.annotationData;
        mutexUnlock();
        return annotationList;
    }

    void MiRNA::getCount(int& n, int& m, int& k, int ann, vector<pair<double, int> >& pvalue)
    {
        m = k = 0;
        n = -1;
        for (int i = 0; i < pvalue.size(); i++) {
            if ( pvalue[i].first >= PERCENT && n < 0 ) n = (int)pvalue.size()-i;
            if ( pvalue[i].second == ann ) {
                if ( n >= 0 ) k++;
                m++;
            }
        }
        return;
    }

    void MiRNA::outputPercent(int num, map<int, int>& index, Matrix& fulldata, vector<vector<PvalData> >& miRNAdata, NMatrix& annotationList)
    {
        for (map<int,int>::iterator it = index.begin(); it != index.end(); it++)
        {
            int count = it->second-1, point = it->first;
            if ( count < 0 ) continue;
            sort( fulldata[count].begin(), fulldata[count].end(), std::greater<double>());
            mutexLock();
            ofstream ofs(pvaluefile.c_str(), ios::app);
            cout << count << endl;
            for( vector<PvalData>::iterator it2 = miRNAdata[count].begin(); it2 != miRNAdata[count].end(); it2++ ) {
                int rank = getRank(it2->pvalue, fulldata[count]), x = point/100, y = point%100;
                //ofs << it2->seed << "\t" << num << "\t" << x << "\t"<< y << "\t" << annotationList[x][y] << "\t" << it2->pvalue << "\t"
                //<< it2->quant << "\t" <<  (double)rank/(double)(fulldata[count]).size() << "\t"
                //<< X.geneSequence->getRandomSeed(it2->seed, GeneSequence::SEEDLENGTH) << endl;
            }
            ofs.close();
            mutexUnlock();
        }
    }
    void MiRNA::getTreeData(int num)
    {
        map<int, int> indexl, indexu;
        Matrix fulldatal, fulldatau;
        vector<vector<PvalData> > miRNAdatal, miRNAdatau;
        setBackGround(false, num, indexl, fulldatal, miRNAdatal);
        setBackGround(true, num, indexu, fulldatau, miRNAdatau);
        NMatrix annotationList = getThreadAnnotationList(num);
        return;
    }

    void MiRNA::setBackGround(bool up, int num, map<int, int>& index, Matrix& fulldata, vector<vector<PvalData> >& miRNAdata)
    {
        Words v;
        int firstseed = -1, count = 1;
        string filename = RawFilename::getpvalueSortFile(num);
        ifstream ifs(filename.c_str());
        fileRead(ifs, v, "\t");
        while( fileRead(ifs, v, "\t") ) {
            if ( v.size() < 9 ) continue;
            else if ( firstseed < 0 ) firstseed = atoi(v[0].c_str());
            int x = atoi(v[2].c_str()), y = atoi(v[3].c_str());
            if ( atoi(v[0].c_str()) == firstseed ) {
                index.insert(pair<int, int>(x*100+y, count++));
                fulldata.resize(count, vector<double>());
                miRNAdata.resize(count, vector<PvalData>());
            }
            double minp = (up) ? atof(v[6].c_str()) : atof(v[5].c_str());
            int minq = (up) ? atoi(v[8].c_str()) : atoi(v[7].c_str());
            if ( !isequal(1,v[1]) ) fulldata[index[x*100+y]-1].push_back(minp);
            else miRNAdata[index[x*100+y]-1].push_back(PvalData(atoi(v[0].c_str()), minp, minq));
        }
    }

    void MiRNA::writePercentThreadMemory(int num)
    {
        NMatrix annotationList = getThreadAnnotationList(num);
        map<int, int> index;
        Matrix fulldata;
        vector<vector<PvalData> > miRNAdata;
        setBackGround(RawFilename::up, num, index, fulldata, miRNAdata);
        outputPercent(num, index, fulldata, miRNAdata, annotationList);
        return;
    }

    void MiRNA::writePercent()
    {
        cout << "write Percent!" << endl;
        ofstream ofs(RawFilename::pvaluefile.c_str());
        ofs << "　seed\tz\tx\ty\tann\tminpvalue\tquantile\tpercent\tsequence" << endl;
        ofs.close();
        vector<ThreadArgMiRNA> thread_arg;
        for (int i = 20; i < 21; i++)
            thread_arg.push_back(ThreadArgMiRNA(i, this));
        thread(1, thread_arg, threadMiRNA);
        return;
    }

     static void *threadMiRNA(void *Arg)
    {
        ThreadArgMiRNA* arg = (ThreadArgMiRNA*)Arg;
        //arg->miRNA->writePercentThread(arg->threadnum);
        arg->miRNA->writePercentThreadMemory(arg->threadnum);
        return NULL;
    }

    template <class Argument>
    void MiRNA::thread(int threadMax, vector<Argument>& thread_arg, void*(*function)(void*))
    {
        pthread_mutex_t *tempmutex = new pthread_mutex_t;
        pthread_mutex_init(tempmutex, NULL);
        mutex = new class Mutex(tempmutex);
        pthread_t threadid[threadMax];
        for (int i = 0; i < threadMax; i++)
            if ( pthread_create(&(threadid[i]), NULL, function, (void*)(&thread_arg[i])) != 0 )
                cout << "error thread" << i << endl;
        for (int i = 0; i < threadMax; i++)
        {
            void *ret = NULL;
            if ( pthread_join(threadid[i], &ret) )
                cout << "error thread end" << i << endl;
        }
        delete mutex;
        mutex = NULL;
        delete tempmutex;
        return;
    }

    void MiRNA::outputAnnRate(int seed, Node *top, vector<int> &numList, vector<int> &allList)
    {
        for (vector<Node*>::iterator it = top->children.begin(); it != top->children.end(); it++)
            outputAnnRate(seed, *it, numList, allList);
        if (top->parent != NULL && allList[top->id] > 0) {
            ofstream ofs(RawFilename::organfile.c_str(), ios::app);
            ofs << seed << "\t" << top->id << "\t" << numList[top->id] << "\t" << allList[top->id] << endl;
            ofs.close();
            numList[top->parent->id] += numList[top->id];
            allList[top->parent->id] += allList[top->id];
            numList[top->id] = allList[top->id] = 0;
        }
    }

    void MiRNA::calcHGD(int seed, vector<PvalPointData>& pList, Node* top)
    {
        vector<int> allList = vector<int>(IDMAX, 0), numList = vector<int>(IDMAX,0);
        sort(pList.begin(), pList.end(), sort_pval);
        int N = (int)pList.size(), n = -1;
        for (int i = 0; i < N; i++) {
            if (n < 0 && NULLHYP/(double)N*(double)i < exp(-pList[i].pvalue)) /*FDR over*/
                n = i;
            if (pList[i].shf == 0) {
                if (n < 0) numList[pList[i].ann]++;
                allList[pList[i].ann]++;
            }
        }
        if (n != 0) {
            if (n < 0) {
                n = N; /* all OK */
                cout << seed << ": FDR 0.05>" << exp(-pList[n-1].pvalue) << " allOK " << n << "/" << N << endl;
            } else
                cout << seed << ": FDR 0.05>" << exp(-pList[n-1].pvalue) << " " << exp(-pList[n].pvalue) << " " << n << "/" << N << endl;
        }
        outputAnnRate(seed, top, numList, allList);
        return;
    }

    void MiRNA::parsePval(int seed, vector<string>& v, vector<PvalPointData>& pList)
    {
        int x = atoi(v[0].c_str()), y = atoi(v[1].c_str()),
            z = atoi(v[2].c_str()), ann = atoi(v[3].c_str());
        for (int i = 4; i < (int)v.size(); i++) {
            double pvalue = atof(v[i].c_str());
            pList.push_back(PvalPointData(seed, i-4, ann, x, y, z, pvalue));
        }
    }

    void MiRNA::readCompressData(int SHFMAX)
    {
        ofstream ofs(RawFilename::organfile.c_str());
        ofs << "seed\tann\tregulationP\tallP " << endl;
        ofs.close();
        AnnotationTree annotation(0);
        annotation.setTree();
        for (int i = 0; i < RawFilename::DIV; i++) {
            string filename = RawFilename::getSeedFile(i);
            ifstream ifs(filename.c_str());
            Words v;
            int seed = 0;
            vector<PvalPointData> pList;
            while( fileRead(ifs, v, "\t") ) {
                if (v.size() == 1) {
                    if (pList.size() > 0) {
                        calcHGD(seed, pList,annotation.first);
                        pList.clear();
                    }
                    seed = atoi(v[0].substr(1).c_str());
                    cout << seed << endl;
                } else if (v.size() == SHFMAX+4)
                    parsePval(seed, v, pList);
            }
            if (pList.size() > 0)
                calcHGD(seed, pList, annotation.first);
        }
        annotation.freeTree();
        return;
    }


    void MiRNA::compressData(int maxthread, int start, int end, int shfmax, int xmax, int ymax)
    {
        Words v;
        vector<PvalList> seedList;
        for (int i = 0; i < MAXSEED; i++)
            seedList.push_back(PvalList(i, shfmax));
        vector<List> annotationData;
        for (int i = 0; i < start; i++)
            annotationData.push_back(List());
        RawFilename::initSeedFile();
        for (int z = start; z < end; z++) {
            class AnnotationTree tree(z);
            annotationData.push_back(tree.getAnnotationNum());
        }
        for (int i = 0; i < maxthread; i++) {
            cout << "--" << i << ":" << std::flush;
            for (int z = start; z < end; z++) {
                cout << z << " " << std::flush;
                string filename = getpvalueFile(z, i);
                ifstream ifs(filename.c_str());
                fileRead(ifs, v, "\t");
                while(fileRead(ifs, v, "\t")) {
                    if (v.size() < 5) continue;
                    seedList[atoi(v[0].c_str())].push(v, z);
                }
            }
            for (vector<PvalList>::iterator it = seedList.begin(); it != seedList.end(); it++) {
                if (it->target) {
                    it->write(annotationData);
                    it->clear();
                }
            }
            cout << endl;
        }
        return;
    }


}
