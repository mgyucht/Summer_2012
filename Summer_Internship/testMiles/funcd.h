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

bool isiMax, isjMax, isjMin;
int iMax, jMax;

struct Funcd {
	
    // The () operator (function call) is now used to calculate the energy of
    // a network. This operator must contain the energy function for the 
    // network.
    
	Doub operator() ( Node ***x, int netSize) {
        
		double funcvalue = 0;
        
        iMax = netSize - 1;
        jMax = netSize - 1;
        
        for(int i=0; i <= iMax; i++) {
            for (int j = 0; j <= jMax; j++) {
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
                
                Node **nPtr = new Node*[4];
                isiMax = i == iMax;
                isjMax = j == jMax;
                isjMin = j == 0;
                
                nPtr[0] = x[i][j];
                
                // PBC is NOT WORKING... :P
                
                int i1 = isiMax ? 0 : i;
                int j1 = isjMax ? 0 : j;
                int j2 = isjMin ? jMax : j;
                
                nPtr[1] = x[i][j1];
                nPtr[2] = x[i1][j];
                nPtr[3] = x[i1][j2];
                
                for (int k = 1; k < 4; k++) {
                    
                    funcvalue += 0.5 * nPtr[k]->sprstiff[k] / nPtr[k]->restlen[k] 
                                    * euclDistSqrd(nPtr[k]->position, nPtr[0]->position, k);
                     
                    printf("%d, %d, %.3f, %.3f, %6.3f | ", i, j, nPtr[k]->position[0], nPtr[0]->position[0], 
                            euclDistSqrd(nPtr[k]->position, nPtr[0]->position, k));
                }
                
                printf("\n");
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
    
    double euclDistSqrd(double *vec1, double *vec2, int k) {
        
        bool xShiftR = false, xShiftL = false, yShiftU = false;
        
        switch (k) {
            case 1: 
                if (isjMax) xShiftR = true;
                break;
            case 2: 
                if (isiMax) yShiftU = true;
                break;
            case 3:
                if (isiMax) yShiftU = true;
                if (isjMax) xShiftL = true;
                break;
        }
        
        double xadj = xShiftR ? jMax * RESTLEN : 0 + xShiftL ? -jMax * RESTLEN : 0;
        double yadj = yShiftU ? iMax * RESTLEN : 0;
        double x = vec2[0] - vec1[0];
        double y = vec2[1] - vec1[1];
        
        return (x - xadj) * (x - xadj) + (y - yadj) * (y - yadj);
    }
};

#endif /*FUNCD_H_*/
