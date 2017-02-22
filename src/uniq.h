//
//  diff.h
//  searchNoise
//
//  Created by carushi on 12/03/29.
//  Copyright 2012年 Undergraduate Student, Department of Science, University of Tokyo. All rights reserved.
//

#ifndef searchNoise_uniq_hpp
#define searchNoise_uniq_hpp

//#include <iostream>
#include <fstream>
#include <string>
#include <vector>
//#include <cstdlib>
//#include <algorithm>

#define LENGTH 50

namespace ushuffle
{
    using std::string;
    using std::sort;
    using std::vector;
    using std::ifstream;
    using std::ofstream;
    using std::endl;
    class Unique
    {
    private:
    public:

        Unique() {}
        ~Unique() {}

    bool uniq(const char* filename, const char* backupfile)
    {
        string str;
        ifstream ifs(filename);
        ofstream ofs(backupfile);
        vector<int> data;
        while ( getline(ifs, str) && str.length() > 0 ) {
            ofs << str << endl;
            data.push_back(atoi(str.c_str()));
        }
        ofs.close(); ifs.close();
        sort(data.begin(), data.end());
        ofs.open(filename);
        int preview = -1;
        for (vector<int>::iterator it = data.begin(); it != data.end(); it++) {
            if ( *it != preview ) ofs << *it << endl;
            preview = *it;
        }
        return true;
    }
    bool uniqFasta(const char* filename, const char* backupfile)
    {
        string str;
        ifstream ifs(filename);
        ofstream ofs(backupfile);
        vector<string> data;
        while ( getline(ifs, str) )
            if ( str.length() > 0 ) {
                ofs << str << endl;
                if ( str[0] == '>' ) data.push_back(str);
            }
        ofs.close(); ifs.close();
        ifs.open(backupfile);
        ofs.open(filename);
        while ( getline(ifs, str) ) {
            if ( str[0] != '>' ) continue;
            vector<string>::iterator it = find(data.begin(), data.end(), str);
            if ( it != data.end() ) {
                ofs << str << endl;
                getline(ifs, str);
                ofs << str << endl;
                data.erase(it);
            }
        }
        return true;
    }
    };

}

#endif

