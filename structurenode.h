//
//  structurenode.h
//  CorrelationClustering
//
//  Created by carushi on 11/11/04.
//  Copyright 2011年 Undergraduate Student, Department of Science, University of Tokyo. All rights reserved.
//

#ifndef CorrelationClustering_structurenode_h
#define CorrelationClustering_structurenode_h

#include <iostream>
#include <string>
#include <vector>
#include <utility>

namespace mousebrain
{
    using std::string;
    using std::vector;
    using std::cout;
    using std::endl;

    class Node
    {
    private:
        int count;
        double r, g, b;

    public:
        int id;
        string _structureName, _abbreviation;
        Node(int Id, double R, double G, double B, string structureName, string abbreviation) :
        id(Id), r(R), g(G), b(B), _structureName(structureName), _abbreviation(abbreviation)
        { count = 0; parent = NULL; children.clear(); }
        ~Node() {}
        Node* parent;
        vector<Node*> children;
        void setParentNode(Node* Parent) { parent = Parent; }
        void freeNode();
        bool upperCount(int Count) { if (count > Count ) return true; else return false; }
        bool isNode(string candidate) { if ( candidate == _abbreviation ) return true; else return false; }
        void printTree(int);
        bool makeChildren(int Id, double R, double G, double B, string structureName, string abbreviation, string parentName);
        bool countNode(int Id);
    };
}


#endif
