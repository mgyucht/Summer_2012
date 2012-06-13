#ifndef _NODE_H_
#define _NODE_H_

// node.h
// ------
//
// node.h provides the data structure for the node segment in the actin network
// simulation. 

#include <iostream>
#include <vector>

using namespace std;

struct CNode {
    
    // dimensions contains the number of dimensions of the system (1 = linear, 
    // 2 = planar, 3 = spatial)
    int dimensions;
    
    // position contains the position of the node segment. sprlen contains
    // the distance between the node segment and a neighboring node segment.
    // sprstiff contains the stiffness modulus of the spring. restlen is the
    // rest length of the spring;
    vector<double> position;
    vector<double> sprlen;
    vector<double> sprstiff;
    double restlen;
    
    // energy contains the current energy of the node segment.
    double energy;
    
    // Cnode(double dimensions, vector<double> position, vector<double> sprlen,
    // vector<double> sprstiff) is the initializer for the node structure. Each
    // of the variables passed is copied into the corresponding variable for
    // the node segment.
    Cnode(int *dim, vector<double> *pos, vector<double> *len,
          vector<double> *stiff, double *rlen) {
        dimensions = *dim;
        position = *pos;
        sprlen = *len;
        sprstiff = *stiff;
        restlen = *rlen;
    }
    
    double calcEnergy() {
        double en = 0;
        
        for (int i = 0; i < (dimensions * (dimensions + 1)) / 2; i++) {
            en += sprstiff[i] / 2 * (sprlen[i] - restlen) * (sprlen[i] - restlen) / restlen;
        }
    }
};

#endif /* _NODE_H_ */
