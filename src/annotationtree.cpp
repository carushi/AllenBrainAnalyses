//
//  annotationtree.cpp
//  CorrelationClustering
//
//  Created by carushi on 11/11/04.
//  Copyright 2011年 Undergraduate Student, Department of Science, University of Tokyo. All rights reserved.
//

#include "annotationtree.h"

namespace mousebrain
{

    void AnnotationTree::makeTree(Node *top)
    {
        ifstream ifs(RawFilename::structureFile.c_str());
        Words v;
        fileRead(ifs, v, ",");
        while( fileRead(ifs, v, ",") ) {
        if ( (int)v.size() < 8 || !top->makeChildren(atoi(v[INFOID].c_str()), atof(v[RED].c_str()), atof(v[BLUE].c_str()), atof(v[GREEN].c_str()), v[STRUCTURE], v[ABBREVIATION], v[PARENT]) )
            cout << "treeerror!" << v[0] << endl;
        }
        ifs.close();
        return;
    }
    void AnnotationTree::setTree()
    {
        first = new class Node(0, 0.0, 0.0, 0.0, "NULL", "NULL");
        first->setParentNode(NULL);
        makeTree(first);
    }

    int AnnotationTree::searchParentId(Node *top, string &name)
    {
        if ( top->_abbreviation == name ) return top->id;
        else
            for (vector<Node*>::iterator it = top->children.begin(); it != top->children.end(); it++)
            {
                int id =  searchParentId(*it, name);
                if ( id >= 0 ) return id;
            }
        return -1;
    }
    int AnnotationTree::getParentId(string &name)
    {
        if ( first == NULL ) setTree();
        if ( name.length() > 0 )
            return searchParentId(first, name);
        return 0;
    }

    void AnnotationTree::cutTree(Node *temp, int limit, vector<int> &realId )
    {
        if ( temp->id != 0 && temp->parent->id != 0 && !temp->upperCount(limit) )
            realId[temp->id] = realId[temp->parent->id];
        for (vector<Node*>::iterator it = temp->children.begin(); it != temp->children.end(); it++)
            cutTree((*it), limit, realId);
        return;
    }
    void AnnotationTree::countAnnotationId(Node *top, List &annotationList)
    {
        ifstream ifs(RawFilename::_annotationDataFileName.c_str());
        Words v;
        fileRead(ifs, v, ","); fileRead(ifs, v, ",");
        while( fileRead(ifs, v, ",") )
        {
            if ( (int)v.size() < 4 ) { /*cout << "error!" << endl;*/ continue; }
            if ( atoi(v[2].c_str()) == z ) {
                int value = atoi(v[3].c_str());
                if ( value < 0 || value >= IDMAX ) continue;
                top->countNode(value);
                annotationList[atoi(v[0].c_str())][atoi(v[1].c_str())] = value;
            }
            else if (atoi(v[2].c_str()) > z ) break;
        }
        ifs.close();
        return;
    }
    void AnnotationTree::simpleCount(List &annotationList)
    {
        ifstream ifs(RawFilename::_annotationDataFileName.c_str());
        Words v;
        fileRead(ifs, v, ","); fileRead(ifs, v, ",");
        while( fileRead(ifs, v, ",") )
        {
            if ( (int)v.size() < 4 ) { /*cout << "error!" << endl;*/ continue; }
            if ( atoi(v[2].c_str()) == z ) {
                int value = atoi(v[3].c_str());
                if ( value < 0 || value >= IDMAX ) continue;
                annotationList[atoi(v[0].c_str())][atoi(v[1].c_str())] = value;
            }
            else if (atoi(v[2].c_str()) > z ) break;
        }
        ifs.close();
        return;
    }

    void AnnotationTree::setParentId(Node *node, vector<int>& pId)
    {
        if (node->children.size() == 0) return;
        for (vector<Node*>::iterator it = node->children.begin(); it != node->children.end(); it++) {
            setParentId(*it, pId);
            pId[(*it)->id] = node->id;
        }
    }

    vector<int> AnnotationTree::getParentIdList()
    {
        setTree();
        vector<int> pId(IDMAX);
        for (int i = 0; i < IDMAX; i++)
            pId[i] = i;
        setParentId(first,pId);
        first->freeNode();
        return pId;
    }

    vector<int> AnnotationTree::getRealId(Node *top, List &annotationList)
    {
        vector<int> realId(IDMAX);
        for (int i = 0; i < IDMAX; i++) realId[i] = i;
        countAnnotationId(top, annotationList);
#ifdef CUT
        cutTree(top, LIMIT, realId);
#endif
        return realId;
    }

    void AnnotationTree::calcRate(Node *top)
    {
        Words v;
        List annotationList(XNUM, vector<int> (YNUM, 0)), clusterElement(ALMOSTCLUSTER, vector<int> (IDMAX, 0));
        vector<int> realId = getRealId(top, annotationList);
        ifstream ifs(RawFilename::_clusteringDataFileName.c_str());
        fileRead(ifs, v, " \"(),");
        while( fileRead(ifs, v, " \"()," ) )
        {
            if ( v.size() < 3 ) continue;
            if ( atoi(v[2].c_str()) > (int)clusterElement.size() ) clusterElement.resize(atoi(v[2].c_str()), vector<int> (IDMAX, 0));
            clusterElement[atoi(v[2].c_str())-1][realId[annotationList[atoi(v[0].c_str())][atoi(v[1].c_str())]]]++;
        }
        ifs.close();
        ofstream ofs(RawFilename::_outputFilename.c_str());
        for (int i = 0; i < (int)clusterElement.size(); i++) {
            ofs << i << "cluster " << std::endl;
            for (int j = 0; j < IDMAX; j++)
                if ( clusterElement[i][j] != 0 ) ofs << j+1 << "infoid ";
            ofs << endl;
            for (int j = 0; j < IDMAX; j++)
                if ( clusterElement[i][j] != 0 ) ofs << clusterElement[i][j] << " ";
            ofs << endl;
        }
        ofs.close();
        return;
    }

    void AnnotationTree::readAnnotation()
    {
        if ( z < 0 || z > 57 ) { cout << "z error!" << endl; return; }
        class Node *top = new class Node(0, 0.0, 0.0, 0.0, "NULL", "NULL");
        top->setParentNode(NULL);
        //cout << "maketree start!" << endl;
        makeTree(top);

        calcRate(top);
        top->freeNode();
        delete top;
        return;
    }

    List AnnotationTree::getAnnotationNum()
    {
        List annotationList(XNUM, vector<int> (YNUM, 0));
        simpleCount(annotationList);
        /*
        class Node *top = new class Node(0, 0.0, 0.0, 0.0, "NULL", "NULL");
        if ( z >= 0 && z < 58 ) {
            top->setParentNode(NULL);
            //cout << "maketree start!" << endl;
            makeTree(top);
            //top->printTree(0);
            vector<int> realId = getRealId(top, annotationList);
            for (int i = 0; i < (int)realId.size(); i++)
                cout << "Id " << i << " " << realId[i] << " ";
            cout << endl;
            for (int i = 0; i < XNUM; i++)
                for (int j = 0; j < YNUM; j++)
                    annotationList[i][j] = realId[annotationList[i][j]];

        }
        top->freeNode();
         */

        return annotationList;
    }
}
