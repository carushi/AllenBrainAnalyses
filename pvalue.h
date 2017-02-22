//
//  pvalue.h
//  searchNoise
//
//  Created by carushi on 12/06/11.
//  Copyright 2012年 Undergraduate Student, Department of Science, University of Tokyo. All rights reserved.
//

#ifndef searchNoise_pvalue_h
#define searchNoise_pvalue_h

//#include <string>
//#include <vector>
#include <map>
#include <fstream>
#include <algorithm>
#include "getstring.h"


namespace mousebrain {
    using std::map;
    using std::string;
    using std::vector;
    using std::ofstream;
    using std::max;
    using std::endl;
    using std::ios;

    typedef vector<vector<int> > List;

    class PvalData
    {
    public:
        PvalData(int Seed, double Pvalue, int Quant) : seed(Seed), pvalue(Pvalue), quant(Quant){}
        ~PvalData(){};
        int seed, quant;
        double pvalue;
    };

    class PvalPointData
    {
    public:
        PvalPointData(int Seed, int Shf, int Ann, int X, int Y, int Z, double Pvalue) : seed(Seed), shf(Shf), ann(Ann), x(X), y(Y), z(Z), pvalue(Pvalue) {}
        ~PvalPointData(){}
        int seed, shf, ann, x, y, z;
        double pvalue;
        bool operator<(const PvalPointData& rhs) {
            return (this->pvalue < rhs.pvalue);
        }
        bool operator>(PvalPointData& rhs) {
            return (this->pvalue > rhs.pvalue);
        }
        bool operator==(const PvalPointData& rhs) {
            return (this->pvalue == rhs.pvalue);
        }
    };


    class PvalList
    {
    private:
        int seed, shfmax;
        map<int, int> _xyzlist;
        vector<vector<double> > _dlist;
    public:
        PvalList(int Seed, int Shfmax) : seed(Seed), shfmax(Shfmax) { target = false; }
        ~PvalList() {}
        bool target;
        int num(int x, int y, int z) { return x+y*70+z*70*70; }

        void set(int key, int &x, int &y, int &z, int &value)
        {
            x = key%70;
            y = (key%4900)/70;
            z = key/4900;
            value = _xyzlist[key];
        }

        void push(vector<string>& v, int z)
        {
            int x = atoi(v[1].c_str()), y = atoi(v[2].c_str()), key = num(x,y,z), shf = atoi(v[5].c_str());
            double pvalue = max(atof(v[3].c_str()), atof(v[4].c_str()));
            if (shf == 1) {
                _xyzlist.insert(map<int, int>::value_type(key, (int)_xyzlist.size()));
                _dlist.push_back(vector<double>(shfmax, -1));
            }
            _dlist[_xyzlist[key]][shf] = pvalue;
            target = true;
        }

        void write(vector<List> &annotationData)
        {
            string filename = RawFilename::getSeedFile(this->seed);
            ofstream ofs(filename.c_str(), ios::app);
            ofs << ">" << seed << endl;
            for (map<int, int>::iterator it = _xyzlist.begin(); it != _xyzlist.end(); it++)
            {
                int x, y, z, value;
                set(it->first, x, y, z, value);
                ofs << x << "\t" << y << "\t" << z << "\t" << annotationData[z][x][y];
                for (vector<double>::iterator it2 = _dlist[value].begin(); it2 != _dlist[value].end(); it2++)
                    ofs << "\t" << *it2;
                ofs << endl;
            }
        }
        void clear()
        {
            _xyzlist.clear();
            _dlist.clear();
            target = false;
        }
    };


}

#endif
