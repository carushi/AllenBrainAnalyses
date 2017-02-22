//
//  genesequence.cpp
//  searchNoise
//
//  Created by carushi on 11/11/21.
//  Copyright 2011年 Undergraduate Student, Department of Science, University of Tokyo. All rights reserved.
//

#include "genesequence.h"


namespace mousebrain
{
    const bool GeneSequence::writemiRNA = false, GeneSequence::debug = false;
    const int DNASeq::ELEMENTNUM = 16, DNASeq::SEEDLENGTH = 7, DNASeq::THREAD = 8,
    DNASeq::MAXMIRNA = 720*2, DNASeq::MAXSHFSEQ = 1000;

    void callSystem(const char* demand, const char* sourcefile, const char* outputfile)
    {
        ostringstream str;
        str << demand << " " << sourcefile << " > " << outputfile;
        system(str.str().c_str());

    }

    void callSystem(const char* demand, const char* sourcefile)
    {
        ostringstream str;
        str << demand << " " << sourcefile;
        system(str.str().c_str());

    }

    string DNASeq::getComplementarySequence(string &seq)
    {
        string str = "";
        reverse(seq.begin(), seq.end());
        for (int i = 0; i < seq.length(); i++)
            str += basePair(seq[i]);
        return str;
    }

    string GeneSequence::getRandomSeed(int count, int length)
    {
        string str;
        int num = MAXSEED/4;
        for (int i = 0; i < length; i++, num /= 4 )
            str.push_back(getChar((count/num)%4));
        return str;
    }

    int GeneSequence::getCount(string& str)
    {
        int sum = 0;
        for (int i = 0; i < str.length(); i++)
            sum = sum*4+encodeChar(str[i]);
        return sum;
    }

    void DNASeq::printNarray(int isid, string& str, int *Narray)
    {
        cout << isid << "Narray " << str << endl;
        cout << isid << "Narray ";
        for (int i = 0; i < THREEBASE; i++)
            cout << Narray[i] << " ";
        cout << endl;
        return;
    }

    string DNASeq::getSeedSequence(ifstream &ifs, int strand, string &name)
    {
        int start = 1, end = 1, state = strand;
        string str = "";
        Words v;
        while ( fileRead(ifs, v, ". \t") )
        {
            if ( v.size() > 0 && istarget(v[0], "//") ) return str;
            if ( v.size() == 4 && istarget(v[0], "FT") && istarget(v[1], "miRNA") ) {
                start = atoi(v[2].c_str());
                end = atoi(v[3].c_str());
                state--;
                if ( state < 0 ) break;
            }
        }
        fileRead(ifs, v, "\""); fileRead(ifs, v, "\"");
        if ( v.size() == 2 ) name = v[1];
        while ( fileRead(ifs, v, ". \t") ) {
            if ( v.size() > 4 && istarget(v[0], "SQ") ) {
                while ( fileRead(ifs, v, " \t") && v.size() > 0 && !istarget(v[0],"//") )
                    for (int i = 0; i < (int)v.size()-1; i++)
                        str += v[i];
                if ( end-start+1 > 0 ) str = str.substr(start-1, end-start+1);
                break;
            }
        }
        str = str.substr(1, SEEDLENGTH);
        str = getComplementarySequence(str);
        return str;
    }

    string DNASeq::getSequence(Words &v)
    {
        ExonList exonSE;
        class Sequence seq;
        Words exonS = myParse(v[9], ","), exonE = myParse(v[10], ",");
        for (int i = 0; i < (int)exonS.size() && i < (int)exonE.size(); i++)
            exonSE.push_back(pair<ulong, ulong>(strtoul(exonS[i].c_str(), NULL, 10)+1, strtoul(exonE[i].c_str(), NULL, 10)));
        string str = "";
        if ( v[6] == v[7] ) return str; /* non coding! */
        else if ( v[3] == "+" ) {
            str = seq.getRawSequence(RawFilename::chrdir, v[2], strtoul(v[7].c_str(), NULL, 10)+1, strtoul(v[5].c_str(), NULL, 10), exonSE);

        } else if ( v[3] == "-" ) {
            str = seq.getRawSequence(RawFilename::chrdir, v[2], strtoul(v[4].c_str(), NULL, 10)+1, strtoul(v[6].c_str(), NULL, 10), exonSE);
            str = getComplementarySequence(str);
        }
        return str;
    }


    void GeneSequence::writeSequence(int fileNum, string &sequence, const char* filename)
    {
        ofstream ofs(filename, ios::app);
        ofs << ">" << fileNum << endl;
        ofs << sequence << endl;
        ofs.close();
    }

    void GeneSequence::writeFileNum(int fileNum, const char* filename)
    {
        ofstream ofs(filename, ios::app);
        ofs << fileNum << endl;
        ofs.close();
    }

     void GeneSequence::writeFilename(vector<int>& filenameList)
    {
        Words v;
        ifstream ifs(RawFilename::filename1.c_str());
        ofstream ofs(RawFilename::filename5.c_str());
        fileRead( ifs, v, "\t" );
        while( fileReadWithoutSkip(ifs, v, "\t") )
        {
            if ( v.size() < 7 || v[4] != "sagittal" ) continue;
            vector<int>::iterator it = find(filenameList.begin(), filenameList.end(), atoi(v[3].c_str()));
            if ( it != filenameList.end() )
                ofs << *it << "\t" << v[6] << endl;
        }
        return;
    }

    void GeneSequence::writeBitSE(int isid, vector<bitset<MAXSEED> >& bsList)
    {
        mutexLock();
        for (int i = 0; i < MAXSHFSEQ+1; i++) {
            writeFileNumBinary(isid, i);
            ofstream ofs(RawFilename::numToSEBF(i).c_str(), ios::binary|ios::app);
            ofs.write((char*)&(bsList[i]), sizeof(bitset<MAXSEED>));
            ofs.close();
        }
        mutexUnlock();

    }
    void GeneSequence::readSeqFromTFile(int threadnum, const char* filename, IDLIST* filenameList)
    {
        ifstream ifs(filename);
        string str;
        int count = 0;
        while (getline(ifs, str)) {
            int isid = getisid(str);
            if ( isid == 0 ) continue;
            bool writeData = (find(filenameList->begin(), filenameList->end(), isid) != filenameList->end()) ? true : false;
            vector<bitset<MAXSEED> > bsList;
            for (int i = 0; i < MAXSHFSEQ+1 && getline(ifs, str); i++)
                if ( writeData ) {
                    if ( i == 0 ) countThreeBaseNum(str,isid);
                    bsList.push_back(countSeedExist(str, isid));
                }
            if (debug) {
                if (writeData) cout << "* ";
                else cout << "# ";
                cout << isid << endl;
            }
            if (!writeData || (int)bsList.size() != MAXSHFSEQ+1 ) continue;
            writeBitSE(isid, bsList);
            count++;
        }
        cout << "--- thread" << threadnum << " writecount " << count << endl;
    }

    bool GeneSequence::stockCheck(const int stock, int threadnum, IDLIST* filenameList)
    {
        if (stock == STOCKMAX) {
            callSystem("/home/username/ushuffle/main", RawFilename::getTFile(threadnum).c_str(), RawFilename::getTFile(threadnum+10).c_str());
            readSeqFromTFile(threadnum, RawFilename::getTFile(threadnum+10).c_str(), filenameList);
            ofstream ofs(RawFilename::getTFile(threadnum).c_str(), ios::trunc); ofs.close();
            return true;
        } else return false;
    }

    void GeneSequence::readSeq(int threadnum, IDLIST* filenameList)
    {
        ifstream ifs(RawFilename::filename2.c_str());
        Words v;
        int stock = 0;
        ofstream ofs(RawFilename::getTFile(threadnum).c_str(), ios::trunc); ofs.close();
        ofs.open(RawFilename::getTFile(threadnum+10).c_str(), ios::trunc); ofs.close();
        ofs.close();
        cout << "--readSeq " << threadnum << endl;
        for (int count = 0; fileRead(ifs, v, " \t"); count++) {
            if ( v.size() >= ELEMENTNUM  && count%THREAD == threadnum-1 ) {
                GeneList::iterator it = imageAndGaId.find(v[1]);
                if ( it != imageAndGaId.end() ) {
                    string sequence = getSequence(v);
                    if ( sequence.length() < MINLENGTH ) continue;
                    changeDNA(sequence);
                    int isid = it->second;
                    mutexLock();
                    writeFileNum(isid, RawFilename::filename7.c_str());
                    writeSequence(isid, sequence, RawFilename::rawsequence.c_str());
                    mutexUnlock();
                    writeSequence(isid, sequence, RawFilename::getTFile(threadnum).c_str());
                    stock++;
                }
                if ( stockCheck(stock, threadnum, filenameList) ) stock = 0;
            }
        }
        stockCheck(STOCKMAX, threadnum, filenameList);

        return;
    }

    void GeneSequence::readFilename()
    {
        ifstream ifs(RawFilename::filename5.c_str());
        Words v;
        imageAndGaId.clear();
        while ( fileRead(ifs, v, "\r\t") )
            if ( v.size() == 2  ) imageAndGaId.insert(pair<string, int> ( v[1], atoi(v[0].c_str()) ));
        cout << "--" << imageAndGaId.size() << " may have sequence" << endl;
        return;
    }

    void GeneSequence::writeMiRNA(int threadnum, int strand)
    {
        ifstream ifs(RawFilename::filename4.c_str());
        Words v;
        for (int count = 0; fileRead(ifs, v, ". \t") && count < MAXMIRNA; )
        {
            if ( v.size() <= 2 ) continue;
            if ( istarget(v[0], "ID") && istarget(v[1].substr(0,3), "mmu") && (++count%THREAD == threadnum-1 )) {
                string name = "", target = getSeedSequence(ifs,strand, name);
                if ( target.length() == SEEDLENGTH ) {
                    if ( name.length() == 0 ) name = v[1];
                    mutexLock();
                    ofstream ofs(RawFilename::miRNAfile.c_str(), ios::app);
                    ofs << name << "\t" << strand << "\t" << target << endl;
                    ofs.close();
                    mutexUnlock();
                }
            }
        }
    }


    void* GeneSequence::singleWrite(MyThreadSeedArg* Arg) /* seed */
    {
        if ( Arg->writeMiRNA ) {
            for (int i = 0; i < 5; i++)
                writeMiRNA(Arg->threadnum, i);
        } else
            readSeq(Arg->threadnum, Arg->filenameList);
        return NULL;
    }

    static void *threadTrans(void *arg)
    {
        MyThreadSeedArg *Arg = (MyThreadSeedArg*)arg;
        GeneSequence *geneSequence = Arg->geneSequence;
        return geneSequence->singleWrite(Arg);
    }

    bitset<MAXSEED> GeneSequence::countSeedExist(string sequence, int fileNum)
    {
        Words seed;
        for (int i = 0; i < MAXSEED; i++)
            seed.push_back(getRandomSeed(i, SEEDLENGTH));
        class SuffixArray sa;
        vector<int> seedNum;
        sa.countSeed(sequence, seed, seedNum);
        bitset<MAXSEED> bs;
        if (seedNum.size() == MAXSEED ) {
            for (int i = 0; i < (int)seedNum.size(); i++)
                if ( seedNum[i] > 0 ) bs.set(i, 1);
        } else
            cout << "seed search error " << fileNum << endl;
        return bs;
    }

    int* GeneSequence::countThreeBaseOfSeed(string &sequence)
    {
        int* Narray = new int[THREEBASE];
        vector<int> count(MAXSEED, 0);
        for (int i = 0, num = 0; i < (int)sequence.length(); i++)
        {
            num = (num)*4+getNum(sequence[i]);
            if ( i >= SEEDLENGTH-1 ) {
                count[num]++;
                num -= getNum(sequence[i-(SEEDLENGTH-1)])*4*4*4*4*4*4;
            }
        }
        for (int i = 0; i < THREEBASE; i++)
            Narray[i] = 0;
        for (int i = 0; i < MAXSEED; i++)
            if ( count[i] > 0 )
                Narray[(i%4)*16+((i/4)%4)*4+(i/16)%4]++;
        return Narray;
    }

    void GeneSequence::countThreeBaseNum(string sequence, int isid)
    {
        class SuffixArray sa;
        int* Narray = sa.countThreeBase(sequence);
        mutexLock();
        ofstream ofs;
        if ( debug && isid < 100 ) {
            ofs.open(getDebugFile(), ios::app);
            ofs << isid << " " << sequence << endl;
            for (int i = 0; i < THREEBASE; i++)
                ofs << Narray[i] << " ";
            ofs << endl;
            ofs.close();
        }
        ofs.open(RawFilename::threebase.c_str(), ios::binary|ios::app);
        ofs.write((char*)&isid, sizeof(int) );
        for (int i = 0; i < THREEBASE; i++)
            ofs.write((char*)&Narray[i], sizeof(int));
        ofs.close();
        geneAnnotate.push_back(isid);
        mutexUnlock();
        delete[] Narray;
    }

    void GeneSequence::throwThread(IDLIST& filenameList)
    {
        vector<MyThreadSeedArg> thread_arg;
        for (int i = 0; i < THREAD; i++) {
            thread_arg.push_back(MyThreadSeedArg(i+1, this));
            thread_arg[i].writeMiRNA = writemiRNA;
            thread_arg[i].setFileList(&filenameList);
        }
        thread::Thread<MyThreadSeedArg> thread;
        thread.variableThread(thread_arg, THREAD, threadTrans);
    }

    void GeneSequence::threadWriteSeedEffect(IDLIST& filenameList)
    {
        if ( RawFilename::local ) return;
        Words v;
        ofstream ofs;
        if (writemiRNA) {
            ofs.open(RawFilename::miRNAfile.c_str(), ios::trunc);
            ofs << "name\tstrand\ttarget" << endl; ofs.close();
            throwThread(filenameList);
        }
        ofs.open(RawFilename::threebase.c_str(), ios::trunc|ios::binary); ofs.close();
        for (int i = 0; i < MAXSHFSEQ+1; i++) {
            ofs.open(RawFilename::numToSEBF(i).c_str(), ios::trunc|ios::binary); ofs.close();
        }
        ofs.open(RawFilename::filename7.c_str(), ios::trunc); ofs.close();
        ofs.open(RawFilename::rawsequence.c_str(), ios::trunc); ofs.close();
        if ( debug ) { ofs.open(getDebugFile(), ios::trunc); ofs.close(); }
        throwThread(filenameList);
        return;
    }

    void GeneSequence::readThreeBaseNum(Data& rankVector, long long int* Narray)
    {
        vector<int> rank;
        for (Data::iterator it = rankVector.begin(); it != rankVector.end(); it++)
            rank.push_back((*it).second);
        ifstream ifs;
        ifs.open(RawFilename::threebase.c_str());
        fill(Narray, Narray+THREEBASE, 0);
        for (int tempid; ifs.read((char*)&tempid, sizeof(int));  ) {
            vector<int>::iterator it = find(rank.begin(), rank.end(), tempid);
            if ( it != rank.end() ) {
                int n;
                for (int i = 0; i < THREEBASE; i++) {
                    ifs.read((char*)&n, sizeof(int));
                    Narray[i] += (long long int)n;
                }
            } else
                for (int i = 0; i < THREEBASE; i++)
                    ifs.read(NULL, sizeof(int));
        }
        return;
    }

    void GeneSequence::addSequence(vector<string>& tList, vector<string>& wList, string seq)
    {
        if (wList.size() == 0) {
            tList.push_back(seq);
        } else {
            for (vector<string>::iterator it = wList.begin(); it != wList.end(); it++)
                tList.push_back(*it+seq);
        }
    }

    vector<int> GeneSequence::getWobble(int count, bool wobble, bool output)
    {
        vector<int> countList;
        if (!RawFilename::wobblepair || !wobble ) {
            countList.push_back(count);
        } else {
            string seed = getRandomSeed(count,SEEDLENGTH);
            vector<string> wobblelist;
            int next;
            for (int i = 0; i < (int)seed.length(); i = next+1) {
                next = getMin((int)seed.find('A',i), (int)seed.find('C',i), (int)seed.length());
                vector<string> tempwoblist;
                addSequence(tempwoblist, wobblelist, seed.substr(i, next-i+1));
                if (next < (int)seed.length()) {
                    string str = (seed[next] == 'A') ? "G" : "T";
                    addSequence(tempwoblist, wobblelist, seed.substr(i, next-i)+str);
                }
                wobblelist = tempwoblist;
            }
            if (output) {
                cout << seed << " " << wobblelist.size() << " : ";
                for (vector<string>::iterator it = wobblelist.begin(); it != wobblelist.end(); it++)
                    cout << *it << " ";
                cout << endl;
            }
            for (vector<string>::iterator it = wobblelist.begin(); it != wobblelist.end(); it++)
                countList.push_back(getCount(*it));
        }
        return countList;
    }

}
