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

const double RESTLEN = 1;

// isiMax, isjMax and isjMin are true whenever i = iMax, j = jMax, or j = jMin
// respectively. Basically, they trigger periodic boundary condition code.

bool isiMax, isjMax, isjMin;

// iMax and jMax refer to the total number of nodes in a single column or in a
// single row, respectively.

int iMax, jMax;

// The Funcd struct.

struct Funcd {
	
    // The () operator (function call) is now used to calculate the energy of
    // a network. This operator must contain the energy function for the 
    // network.
    //
    // Node ***x is a Node two-dimensional array, and the size of that array is
    // netSize x netSize.
    
	Doub operator() ( Node ***x, int netSize, int strain ) {
        
		double funcvalue = 0;
        
        iMax = netSize - 1;
        jMax = netSize - 1;
        
        for(int i = 0; i <= iMax; i++) {
            for (int j = 0; j <= jMax; j++) {
                
                // Make the pointer array Node *nPtr[]. This shall point to the
                // following nodes:
                //  
                //  1. node i, j
                //  2. node i, j + 1 (j1)
                //  3. node i + 1 (i1), j
                //  4. node i + 1, j - 1 (j2)
                //  
                // Make sure to impose the periodic boundary condition here. 
                // This is done by making a temporary node copy of the 
                // corresponding mode (using modulo arithemetic). 
                //
                // Strain is obtained by shifting the bonds at the top of the
                // network over by an integer number of points in the
                // x-direction. When measuring the energy, it only affects the
                // topmost row (i = iMax).
                
                Node **nPtr = new Node*[4];
                isiMax = i == iMax;
                isjMax = j == jMax;
                isjMin = j == 0;
                
                nPtr[0] = x[i][j];
                
                int j1 = isjMax ? 0 : j + 1;
                int i1 = isiMax ? 0 : i + 1;
                int j2 = isjMin ? jMax : j - 1;
                
                nPtr[1] = x[i][j1];
                nPtr[2] = x[i1][j];
                nPtr[3] = x[i1][j2];
                
                for (int k = 1; k < 4; k++) {
                    
                    funcvalue += 0.5 * nPtr[k]->sprstiff[k] / nPtr[k]->restlen[k]
                                    * deltaLSqrd(nPtr, k, strain);
                    
                }
            } // Cycling all columns in a row
		} // Cycling all rows
        
		return funcvalue;
	}
	
    // df(Node ***x, double ***dx, int netSize) is the gradient function for the 
    // network. This determines the gradient function for a point in the
    // network.
    
    void df(Node ***x, double ***dx, int netSize) {
        
        iMax = netSize - 1;
        jMax = netSize - 1;
        
        for (int i = 0; i <= iMax; i++) {
            for (int j = 0; j <= jMax; j++) {
                //
                // TODO:GRADIENT FUNCTION AUGH
                //
            }
        }
    }

private: 
    
    // deltaLSqrd returns the square change in the length of the spring
    // between nodes nPtr[0] and nPtr[k].
    
    double deltaLSqrd(Node **x, int k, int strain) {
        double lijk = euclDistSqrd(x[0]->position, x[k]->position, k, strain);
        double l0 = RESTLEN;
        
        return (lijk - l0) * (lijk - l0);
    }
    
    // euclDistSqrd is a Euclidean distance calculator that works with the
    // pointer nPtr declared in operator(). It takes into account the periodic
    // boundary condition for the system. This is the function that incorporates
    // strain into the network.
    // 
    // Note that we use jMax + 1 and iMax + 1 in these expressions because they
    // are equivalent to the size of the network netSize. 
    
    double euclDistSqrd(double *vec1, double *vec2, int k, int strain) {
        
        double xshift = 0.0, yshift = 0.0;
        
        switch (k) {
            
            // Case where nPtr[1] is in column 0.
            case 1: 
                if (isjMax)
                    
                    xshift += (jMax + 1) * RESTLEN;
                
                break;
                
            // Case where nPtr[2] is in row 0.
            case 2: 
                if (isiMax) {
                
                    xshift += (jMax + strain + 1) * RESTLEN / 2.0;
                    yshift += (iMax + 1) * RESTLEN * sqrt(3.0) / 2.0;
                    
                }
                break;
                
            // Case where nPtr[3] is in row 0 and/or column jMax.
            case 3:
                if (isiMax) {
                
                    xshift += (jMax + strain + 1) * RESTLEN / 2.0;
                    yshift += (iMax + 1)* RESTLEN * sqrt(3.0) / 2.0;
                    
                }
                if (isjMin)
                    
                    xshift -= (jMax + 1) * RESTLEN;

                break;
        }
        
        // xadj and yadj refer to nPtr[k = 1, 2, 3].
        
        double xadj = vec2[0] + xshift; 
        double yadj = vec2[1] + yshift;
        
        // x and y refer to nPtr[0] always.
        
        double x = vec1[0];
        double y = vec1[1];
        
        return (xadj - x) * (xadj - x) + (yadj - y) * (yadj - y);
    }
};

#endif /*FUNCD_H_*/
