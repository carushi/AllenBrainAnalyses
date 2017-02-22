//
//  burrowswheelertransform.cpp
//  searchNoise
//
//  Created by carushi on 11/12/02.
//  Copyright 2011年 Undergraduate Student, Department of Science, University of Tokyo. All rights reserved.
//

#include "burrowswheelertransform.h"

namespace mousebrain
{
    int BWT::getOcc(int c, int limit, Array& Occ)
    {
        int new_limit = 0;
        if ( limit < slen ) {
            new_limit = Occ[limit/OCCSIZE][c-1];
            for (int i = limit-limit%OCCSIZE+1; ( limit%OCCSIZE ) && i <= limit; i++) {
                int j = ( SA[i] ) ? SA[i]-1 : slen-1;
                if ( s[j] == c ) new_limit++;
            }
        }
        return new_limit;
    }

    int BWT::search(int *query, int length, Array &Occ, bool up)
    {
        int limit = ( up ) ? 0 : slen-1;
        for (int i = 0; i < length; i++)
        {
            int c = query[length-1-i];
            if ( up && !limit ) limit = slen+1;
            limit = ( up ) ? bkt[c-1]+getOcc(c, limit-1, Occ) : bkt[c-1]+getOcc(c, limit, Occ)-1;
        }
        return limit;
    }

    void BWT::makeOcc(Array &Occ)
    {
        fill(bkt, bkt+K, 0);
        for (int i = 0; i < slen; i++)
        {
            int c = ( SA[i] == 0 ) ? s[slen-1] : s[SA[i]-1];
            if ( c != 0 )	bkt[c-1]++;
            if ( i%OCCSIZE == 0 )
                for ( int j = 0; j < K; j++ )
                    Occ[i/OCCSIZE][j] = bkt[j];
        }
        for ( int i = K; i > 0; i-- )
        {
            bkt[i-1] = 0;
            for ( int j = i; j > 1; j-- )
                bkt[i-1] += bkt[j-2];
            bkt[i-1]++;
        }
    }

    void BWT::getQuery(string &str, int *query)
    {
        for (int i = 0; i < str.length(); i++)
        {
            if (str[i] == '$') query[i] = 0;
            else if (str[i] == 'C' || str[i] == 'c') query[i] = 2;
            else if (str[i] == 'G' || str[i] == 'g') query[i] = 3;
            else if (str[i] == 'T' || str[i] == 't') query[i] = 4;
            else query[i] = 1;
        }
        return;
    }

    void BWT::printDEBUG(Array &Occ, int u, int l, int *query, int length)
    {
        /*
        for	( int i = 0; i < (slen-1)/OCCSIZE; i++ ) {
            for ( int j = 0; j < K; j++ )
                cout << Occ[i][j] << " ";
        cout << endl;
        }
        cout << "bkt";
        for (int i = 0; i < K; i++)
            cout << bkt[i] << " ";
        for (int i = 0; i < length; i++)
            cout << query[i] << " ";
        cout << endl;
         */
        for (int i = l; i < u; i++)
        {
            for (int j = 0; j < length; j++)
            {
                if (SA[i]+j == slen) break;
                cout << s[SA[i]+j] << " ";
            }
            cout << endl;
        }
        return;
    }

    void BWT::getSeedNum(vector<string>& seed, vector<int>& seedNum)
    {
        if (seed.size() == 0) return;
        int seedlen = (int)seed[0].length(), *query = new int[seedlen];
        Array Occ(( slen-1 )/OCCSIZE+1, vector<int>(K));
        makeOcc(Occ);
        for (vector<string>::iterator it = seed.begin(); it != seed.end(); it++) {
            getQuery(*it, query);
            int l = search(query, seedlen, Occ, true ),
            u = search(query, seedlen, Occ, false );
            if ( l <= u ) seedNum.push_back(u-l+1);
            else seedNum.push_back(0);
        }
       delete[] query;
    }

    int BWT::BWTransform(string &str)
    {
        int *query = new int[str.length()];
        getQuery(str, query);
        Array Occ(( slen-1 )/OCCSIZE+1, vector<int>(K));
        makeOcc(Occ);

        int l = search(query, (int)str.length(), Occ, true ),
            u = search(query, (int)str.length(), Occ, false );
        delete[] query;
        if ( l <= u ) return u-l+1;
        else return 0;
    }

    void BWT::getQueryThreeBase(int *query, int num)
    {
        query[0] = num % 4+1;
        query[1] = (num / 4) % 4+1;
        query[2] = num / 16+1;
    }
     int* BWT::getThreeBaseCount()
    {
        int *query = new int[3];
        Array Occ(( slen-1 )/OCCSIZE+1, vector<int>(K));
        makeOcc(Occ);

        int* Narray = new int[THREEBASE];
        for (int i = 0; i < THREEBASE; i++)
        {
            getQueryThreeBase(query, i);
            int l = search(query, 3, Occ, true),
            u = search(query, 3, Occ, false);
            //cout << query[0] << query[1] << query[2]  << " " << l << " " << u << endl;
            //printDEBUG(Occ, u, l, query, 3);
            if (RawFilename::mask) {
            vector<int> startSite;
            for (int j = l; j <= u; j++) startSite.push_back(SA[j]);
            sort(startSite.begin(), startSite.end());
            int preview = -1;
            for (vector<int>::iterator it = startSite.begin(); it != startSite.end(); it++)
                if ( preview < 0 || preview+2 < *it ) preview = *it;
                else u--;
            }
            Narray[i] = ( l <= u ) ? u-l+1 : 0;
        }
        delete[] query;
        return Narray;
    }

}
