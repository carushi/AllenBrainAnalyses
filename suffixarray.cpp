//
//  suffixarray.cpp
//  searchNoise
//
//  Created by carushi on 11/11/29.
//  Copyright 2011年 Undergraduate Student, Department of Science, University of Tokyo. All rights reserved.
//

#include "suffixarray.h"
namespace mousebrain
{
    void SuffixArray::printSA(int Length)
    {
        cout << "printSA (" << Length << ")";
        for (int i = 0; i < Length; i++)
            cout << SA[i] <<  " ";
        cout << endl;
    }

    void SuffixArray::prints(int *s, int Length)
    {
        cout << "prints ";
        for (int i = 0; i < Length; i++)
            cout << s[i] << " ";
        cout << endl;
    }

    void SuffixArray::reserveArray(int Length, string& sequence)
    {
        length = Length;
        SA = new int[length];
        fill(SA, SA+length, 0);
        str = new int[length];
        for (int i = 0; i < Length; i++)
            str[i] = encodeChar(sequence[i]);
        return;
    }

    void SuffixArray::freeArray()
    {
        delete[] SA;
        delete[] str;
        return;
    }

    void SuffixArray::getBucket(int *s, int *bkt, int Length, int K, bool start)
    {
        int sum = 0;
        fill(bkt, bkt+K, 0);
        for ( int i = 0; i < Length; i++ )
            bkt[s[i]]++;
        for ( int i = 0; i < K; i++) {
            sum += bkt[i];
            bkt[i] = (start) ? sum : sum-bkt[i];
        }
        return;
    }

    void SuffixArray::inducedSorting(int *s, int *t, int *bkt, int Length, int K, bool ltype)
    {
        getBucket(s, bkt, Length, K, !ltype);
        if ( ltype ) {
            for ( int i = 0; i < Length; i++ ) {
                if ( SA[i]-1 >= 0 && !tget(t, SA[i] - 1, Length) )
                        SA[bkt[s[SA[i]-1]]++] = SA[i]-1;
            }
        } else {
            for ( int i = Length-1; i >= 0; i--)
                if ( SA[i]-1 >= 0 && tget(t, SA[i] - 1, Length) )
                        SA[--bkt[s[SA[i]-1]]] = SA[i]-1;
        }
        return;
    }

    void SuffixArray::tempSetting(int *s, int *bkt, int lms)
    {
        for (int i = lms-1; i >= 0; i--) {
            int temp = SA[i];
            SA[i] = -1;
            SA[--bkt[s[temp]]] = temp;
        }
        return;
    }

    bool SuffixArray::moveLMS(int* s, int *t, int *bkt, int Length, int lms)
    {
        int name = 0, prev = -1;
        for (int i = 0; i < lms; i++) {
            int pos = SA[i];
            bool diff = false;
            for (int d = 0; d < Length; d++) {
                if (prev == -1 || s[pos+d] != s[prev+d] || tget(t, pos+d, Length) != tget(t, prev+d, Length)) {
                    diff = true; break;
                } else if (d > 0 && isLMS(t, pos+d, Length)) break;
            }
            if (diff) { name++; prev = pos; }
            pos = ( pos%2 == 0 ) ? pos/2:(pos-1)/2;
            SA[lms+pos] = name-1;
        }
        for (int i = Length-1, j = Length-1; i >= lms; i-- )
            if ( SA[i] >= 0 )
                    SA[j--] = SA[i];
        if ( name < lms ) {
            SAIS(SA+Length-lms, lms, name);
            return false;
        } else return true;
    }

    void SuffixArray::SAIS(int *s, int Length, int K)
    {
        int *t = new int[Length+1], *bkt = new int[K], lms = 0;
        tset(t, Length-1, 1, Length); tset(t, Length-2, 0, Length);
        for (int i = Length - 3; i >= 0; i--) tset(t, i, ( s[i] < s[i+1] || ( s[i] == s[i+1] && tget(t, i+1, Length) != 0 ) ) ? 1 : 0, Length );
        getBucket(s, bkt, Length, K, true);
        fill(SA, SA+Length, -1);
        for ( int i = 1; i < Length; i++ )
            if ( isLMS(t, i, Length) ) SA[--bkt[s[i]]] = i;
        inducedSorting(s, t, bkt, Length, K, true);
        inducedSorting(s, t, bkt, Length, K, false);
        for (int i = 0; i < Length; i++)
            if ( isLMS(t, SA[i], Length) ) SA[lms++] = SA[i];
        fill(SA+lms, SA+Length, -1);
        if ( moveLMS(s, t, bkt, Length, lms) )
            for ( int i = 0; i < lms; i++ )
                    SA[SA[Length-lms+i]] = i;
        getBucket(s, bkt, Length, K, true);
        for (int i = 1, j = 0; i < Length; i++) {
            if ( isLMS( t, i, Length ) )
				SA[Length-lms+j++] = i;
		}
        for ( int i = 0; i < lms; i++ )
            SA[i] = SA[Length-lms+SA[i]];
        fill(SA+lms, SA+Length-1, -1);
        tempSetting(s, bkt, lms);
        inducedSorting(s, t, bkt, Length, K, true);
        inducedSorting(s, t, bkt, Length, K, false);
        delete[] t;
        delete[] bkt;

        return;
    }

    void SuffixArray::printDebug(string &seed)
    {
        for (int i = 0; i < length; i++)
        {
            cout << i << " ";
            for (int j = 0; j < seed.length() && SA[i]+j < length; j++)
                cout << str[SA[i]+j] << " ";
            cout << endl;
        }
        return;
    }

    int SuffixArray::countSequence(string &sequence, string &seed)
    {
        int num = 0;
        sequence.append("$");
        reserveArray((int)sequence.length(), sequence);
        SAIS(str, length, BASE);
        if ( seed.length() < sequence.length() ) {
            class BWT bwt(str, SA, length, BASE);
            num = bwt.BWTransform(seed);
        }
        freeArray();
        return num;
    }

    int* SuffixArray::countThreeBase(string &sequence)
    {
        sequence.append("$");
        reserveArray((int)sequence.length(), sequence);
        SAIS(str, length, BASE);
        class BWT bwt(str, SA, length, BASE);
        int* Narray = bwt.getThreeBaseCount();
        freeArray();
        return Narray;
    }

    string SuffixArray::countSeed(string &sequence, vector<string> &seed)
    {
        sequence.append("$");
        reserveArray((int)sequence.length(), sequence);
        SAIS(str, length, BASE);
        class BWT bwt(str, SA, length, BASE);
        vector<int> count;
        for (vector<string>::iterator it = seed.begin(); it != seed.end(); it++)
            count.push_back(bwt.BWTransform(*it));
        freeArray();
        return RawFilename::countToStr(count);
    }

    void SuffixArray::countSeed(string& sequence, vector<string>&seed, vector<int>& seedNum)
    {
        sequence.append("$");
        reserveArray((int)sequence.length(), sequence);
        SAIS(str, length, BASE);
        class BWT bwt(str, SA, length, BASE);
        bwt.getSeedNum(seed, seedNum);
        freeArray();
        return;
    }

}
