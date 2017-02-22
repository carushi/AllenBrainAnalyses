//
//  annotationtree.h
//  CorrelationClustering
//
//  Created by carushi on 11/11/04.
//  Copyright 2011年 Undergraduate Student, Department of Science, University of Tokyo. All rights reserved.
//

#ifndef CorrelationClustering_annotationtree_h
#define CorrelationClustering_annotationtree_h

//#include <cstdlib>
#include "structurenode.h"
#include "getstring.h"

#define STRUCTURE (0)
#define ABBREVIATION (1)
#define PARENT (2)
#define RED (3)
#define GREEN (4)
#define BLUE (5)
#define INFOID (6)
#define STRID (7)
#define LIMIT (21)
#define IDMAX (233)
#define ALMOSTCLUSTER (400)
//#define CUT

namespace mousebrain
{
    using std::cout;
    using std::endl;
    using std::ofstream;

    typedef vector<vector<int> > List;

    class AnnotationTree
    {
    private:
        int z;
        void makeTree(Node *);
        void calcRate(Node *);
        void cutTree(Node *, int, vector<int> &);
        void countAnnotationId(Node *, List &);
        void simpleCount(List &);
        void readAnnotation();
        int searchParentId(Node *, string &);
    public:
        AnnotationTree(int Z) : z(Z) {
            first = NULL;
        }
        ~AnnotationTree() {}
        Node *first;
        static const int XNUM = 67, YNUM = 41, ZNUM = 58;
        void setTree();
        void freeTree() { first->freeNode(); }
        vector<int> getRealId(Node *top, List &annotationList);
        List getAnnotationNum();
        void setParentId(Node *, vector<int>&);
        int getParentId(string &);
        vector<int> getParentIdList();

    };
}

#endif
