// funcd.h
// -------
// 
// funcd.h provides the energy function and derivative functions for a network.
// This header file is designed to be used for the conjugate gradient method 
// of minimizing the value of a multidimensional function. 
//
// operator() ( VecNode2D &x )
// ---------------------------
//
// Reports the energy of the network x. 
// 
// df ( VecDoub_I &x, VecDoub_O &deriv)
// ------------------------------------
//
// Calculates the gradient of the network's energy.

#ifndef FUNCD_H_
#define FUNCD_H_

#include "frprmn.h"
#include "node.h"

// length must be even: # DOF's along x will be 2*(dimx-1)
const int length=40;

const double RESTLEN = 1;

struct Funcd {
	
    // The () operator (function call) is now used to calculate the energy of
    // a network. This operator must contain the energy function for the 
    // network.
    
	Doub operator() ( VecNode2D &x ) {
		double funcvalue=0;
        
        for(int i=0; i < x.nrows(); i++) {
            for (int j = 0; j < x.ncols(); j++) {
                // Make the pointer array Node *nPtr[]. This shall point to the
                // following nodes:
                //  
                //  1. node i, j
                //  2. node i, j + 1
                //  3. node i + 1, j
                //  4. node i + 1, j - 1
                //  
                // Make sure to impose the periodic boundary condition here. 
                // This is done by making a temporary node copy of the 
                // corresponding mode (using modulo arithemetic). This should
                // also take into account strain.
            }
		}
        
		return funcvalue;
	}
	
    // df(VecDoub_I &x, VecDoub_O &deriv) is the gradient function for the 
    // network. This determines the gradient function for a point in the
    // network.

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
