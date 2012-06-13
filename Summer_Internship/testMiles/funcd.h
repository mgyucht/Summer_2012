// funcd.h
// -------
// 
// funcd.h provides the energy function and derivative functions for a network.

#ifndef FUNCD_H_
#define FUNCD_H_

#include "frprmn.h"
#include "node.h"

const double MU = 1;
const double RESTLEN = 1;

// length must be even: # DOF's along x will be 2*(dimx-1)
const int length=40;

struct Funcd {
	
    // The () operator (function call) is now used to calculate the energy of
    // a network. This operator must contain the energy function for the 
    // network.
    
	Doub operator() ( VecNode2D &x ) {
		double funcvalue=0;
        
        for(int i=0; i < x.nrows(); i++) {
            for (int j = 0; j < x.ncols(); j++) {
                funcvalue += x[i][j].calcEnergy();
            }
		}
        
		return funcvalue;
	}
	
    // df(VecDoub_I &x, VecDoub_O &deriv) is the gradient function for the 
    // network. This determines 

	void df(const VecDoub_I &x, VecDoub_O &deriv) 
	{
        for (int i=0; i<length; i++) {
            // INSERT GRADIENT FUNCTION HERE
		}
	}
    
    void df(const VecNode2D &x, VecNode2D &deriv) {
        for (int i = 0; i < x.nrows(); i++) {
            for (int j = 0; j < x.ncols(); j++) {
                // GRADIENT FUNCTION HERE
            }
        }
    }
};

#endif /*FUNCD_H_*/
