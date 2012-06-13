#ifndef _NODE_H_
#define _NODE_H_

// node.h
// ------
//
// node.h provides the data structure for the node segment in the actin network
// simulation. 

#include <iostream>
#include <vector>
#include "nr3.h"

using namespace std;

struct Node {
    
    // position contains the position of the node segment. sprlen contains
    // the distance between the node segment and a neighboring node segment.
    // sprstiff contains the stiffness modulus of the spring. restlen is the
    // rest length of the spring. 
    //
    // In these vectors, the indexing is as follows: for position, the first
    // entry is the abcissa, and the second entry is the ordinate. For sprlen
    // and sprstiff for a node at (i, j), the first entry corresponds to the
    // bond with node at (i, j + 1), the second entry with node at (i + 1, j), 
    // and third entry with node at (i + 1, j - 1).
    
    vector<double> position;
    vector<double> sprlen;
    vector<double> sprstiff;
    double restlen;
    
    // energy contains the current energy of the node segment.
    
    double energy;
    
    // Node(vector<double> position, vector<double> sprlen,
    // vector<double> sprstiff) is the initializer for the node
    // structure. Each of the variables passed is copied into the
    // corresponding variable for the node segment.
    
    Node(vector<double> pos, vector<double> len, vector<double> stiff, 
          double rlen) {
        position = pos;
        sprlen = len;
        sprstiff = stiff;
        restlen = rlen;
    }
    
    double calcEnergy() {
        double en = 0;
        
        for (int i = 0; i < sprlen.size(); i++) {
            en += sprstiff[i] / 2 * (sprlen[i] - restlen) * (sprlen[i] - restlen) / restlen;
        }
        
        return en;
    }
};

typedef NRvector<Node> VecNode1D;
typedef NRmatrix<Node> VecNode2D;
typedef NRmat3D<Node> VecNode3D;
typedef VecNode1D VecNode;

#endif /* _NODE_H_ */
