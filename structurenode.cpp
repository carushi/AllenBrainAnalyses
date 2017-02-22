//
//  structurenode.cpp
//  CorrelationClustering
//
//  Created by carushi on 11/11/04.
//  Copyright 2011年 Undergraduate Student, Department of Science, University of Tokyo. All rights reserved.
//

#include "structurenode.h"

namespace mousebrain
{
    void Node::freeNode()
    {
        for (vector<Node*>::iterator it = children.begin(); it != children.end(); it++)
        {
            (*it)->freeNode();
            delete (*it);
        }
        return;
    }

    void Node::printTree(int depth)
    {
        /*
        cout << "id:" << id;
        cout << " " << _abbreviation;
        cout << " children:" << children.size();
        cout << " count:" << count << " ";
        if ( parent != 0 ) cout << "parent:" << (*parent).id << endl;
        else cout << "parent: NULL" << endl;
        for (vector<Node*>::iterator it = children.begin(); it != children.end(); it++)
            (*it)->printTree();
         */
        cout << id << " " << depth << endl;
        for (vector<Node*>::iterator it = children.begin(); it != children.end(); it++)
            (*it)->printTree(depth+1);

        return;
    }

    bool Node::makeChildren(int Id, double R, double G, double B, string structureName, string abbreviation, string parentName)
    {
        if ( parentName == "NULL" || isNode(parentName) ) {
            class Node *new_child = new Node( Id, R, G, B, structureName, abbreviation );
            new_child->setParentNode(this);
            children.push_back(new_child);
            return true;
        }
        for (vector<Node*>::iterator it = children.begin(); it != children.end(); it++)
            if ( (*it)->makeChildren(Id, R, G, B, structureName, abbreviation, parentName) ) return true;
        return false;
    }

    bool Node::countNode(int Id)
    {
        if ( Id == id ) {
            this->count++;
            for (Node* temp = this->parent; temp->id != 0; temp = temp->parent) (*temp).count++;
            return true;
        }
        for (vector<Node*>::iterator it = children.begin(); it != children.end(); it++)
            if ( (*it)->countNode(Id) == true ) return true;
        return false;
    }
}
