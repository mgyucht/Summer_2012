#ifndef _NODE_H_
#define _NODE_H_

// node.h
// ------
//
// node.h provides the data structure for the node segment in the actin network
// simulation. 

#include <iostream>
#include "nr3.h"

using namespace std;

struct Node {
    
    // position contains the position of the node segment. 
    // sprstiff contains the stiffness modulus of the spring. restlen is the
    // rest length of the spring. 
    //
    // In these vectors, the indexing is as follows: for position, the first
    // entry is the abcissa, and the second entry is the ordinate. For restlen and
    // sprstiff for a node at (i, j), the first entry corresponds to the
    // bond with node at (i, j + 1), the second entry with node at (i + 1, j), 
    // and third entry with node at (i + 1, j - 1).
    
    double *position, *sprstiff, restlen;
    
    // Node(vector<double> position, vector<double> sprlen,
    // vector<double> sprstiff) is the initializer for the node
    // structure. Each of the variables passed is copied into the
    // corresponding variable for the node segment.
    
    Node() {
        position = new double[2];
        sprstiff = new double[3];
        restlen = 0;
    }
    
    Node(double *pos, double *stiff, double rlen) : position(pos), 
            sprstiff(stiff), restlen(rlen) {}
    
    Node(Node & a) {
        
        this->position = a.position;
        this->sprstiff = a.sprstiff;
        this->restlen = a.restlen;
        
    }

    ~Node() {
        delete[] position;
        delete[] sprstiff;
    }
    
    // printPos() is designed to test the functionality of the Node struct.
    // It prints the position as a coordinate in the xy-plane. 
    // Sample output:
    // (0, 0.3)
    
    void printPos() {
        
        double xpos = position[0];
        double ypos = position[1];

        cout << setw(5) << "(" << xpos << ", " << ypos << ") ";
    }
    
};

typedef NRvector<Node> VecNode1D;
typedef NRmatrix<Node> VecNode2D;
typedef NRMat3d<Node> VecNode3D;
typedef VecNode1D VecNode;

#endif /* _NODE_H_ */
