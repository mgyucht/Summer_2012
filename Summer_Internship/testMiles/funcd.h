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

#include <vector>
#include <stdio.h>
#include "debug.h"

#define xadj pos[2 * k] + xshift 
#define yadj pos[2 * k + 1] + yshift

using namespace std;

const double RESTLEN = 1;
const double DEL = 1e-11;

// isiMax, isjMax and isjMin are true whenever i = iMax, j = jMax, or j = jMin
// respectively. Basically, they trigger periodic boundary condition code.

bool isiMax, isiMin, isjMax, isjMin;

// iMax and jMax refer to the total number of nodes in a single column or in a
// single row, respectively.

int iMax, jMax;

extern int netSize;
extern double strain;

// The Funcd struct.

struct Funcd {
    
    double *** spr; 
    double *** rlen;
    int i1, i2, j1, j2;
    
    // Constructor. Taking the correct approach.

    Funcd(double *** ssprstiff, double *** rrestlen) : spr(ssprstiff), rlen(rrestlen) {}

    // The () operator (function call) is now used to calculate the energy of
    // a network. This operator must contain the energy function for the
    // network.
    //
    // double *x is the position array. Information in this array is arranged
    // accordingly:
    //
    // pos[0] => Node (0, 0), abcissa
    // pos[1] => Node (0, 0), ordinate
    // pos[2] => Node (1, 0), abcissa
    // pos[3] => Node (1, 0), ordinate
    // ...
    // pos[2 * netSize] => Node (0, 1), abcissa
    // pos[2 * netSize + 1] => Node (0, 1), ordinate
    // pos[2 * netSize + 2] => Node (1, 1), abcissa
    // pos[2 * netSize + 3] => Node (1, 1), ordinate

    double operator() ( double *x ) {
        
        #if DEBUGENERGY
        
            printf("Energy function called. \n");
            
        #endif
        
        double funcvalue = 0;

        iMax = netSize - 1;
        jMax = netSize - 1;

        for(int i = 0; i <= iMax; i++) {
            for (int j = 0; j <= jMax; j++) {

                // Make the array double pos[]. This shall point to the
                // positions of the following nodes:
                //
                // 1. node i, j
                // 2. node i, j + 1 (j1)
                // 3. node i + 1 (i1), j
                // 4. node i + 1, j - 1 (j2)
                //
                // Make sure to impose the periodic boundary condition here.
                // This is done by making a temporary node copy of the
                // corresponding mode (using modulo arithemetic).
                //
                // Strain is obtained by shifting the bonds at the top of the
                // network over by an integer number of points in the
                // x-direction. When measuring the energy, it only affects the
                // topmost row (i = iMax).

                isiMax = i == iMax;
                isjMax = j == jMax;
                isjMin = j == 0;

                j1 = isjMax ? 0 : j + 1;
                i1 = isiMax ? 0 : i + 1;
                j2 = isjMin ? jMax : j - 1;
                
                double pos[8] = {
                    
                    x[(i * netSize + j) * 2], 
                    x[(i * netSize + j) * 2 + 1], 
                    x[(i * netSize + j1) * 2],
                    x[(i * netSize + j1) * 2 + 1],
                    x[(i1 * netSize + j) * 2],
                    x[(i1 * netSize + j) * 2 + 1],
                    x[(i1 * netSize + j2) * 2],
                    x[(i1 * netSize + j2) * 2 + 1]
                        
                };
                
                #if DEBUGENERGY
                
                    printf("(%d, %d):\n", j, i);
                    
                #endif
                
                for (int k = 1; k < 4; k++) {
                    
                    double temp = 0.5 * spr[i][j][k - 1] / rlen[i][j][k - 1] 
                                    * deltaLSqrd(pos, i, j, k);
                    
                    #if DEBUGENERGY
                    
                        printf("    %d: (%.2f, %2f), %3f\n", k, pos[2 * k], 
                                pos[2 * k + 1], temp);
                        
                    #endif
                    
                    funcvalue += temp;
                
                }
                    
            } // Cycling all columns in a row
        } // Cycling all rows
        
        return funcvalue;
    }

    // df(double *x, double *dx) is the gradient function for the
    // network. This determines the gradient function for a point in the
    // network.

    void df(double *x, double *dx) {
        
        #if DEBUGGRADIENT
            
        printf("Gradient function called.\n");
        
        #endif
        
        // dispdx is the Funcd object used to calculate the energy of the
        // original network pointed by x and the network pointed by dispx.
        // dispx is identical to x in all entries except for one, where it
        // differs by DEL declared at the top of this file. The gradient is then
        // calculated by the change in energy associated with this change in the
        // array.

        iMax = netSize - 1;
        jMax = iMax;
        
        
        for (int i = 0; i <= iMax; i++) {
            for (int j = 0; j <= jMax; j++) {
                
                double current, dxtemp, dytemp;
                current = dxtemp = dytemp = 0.0;
                
                isiMax = i == iMax;
                isjMax = j == jMax;
                isiMin = i == 0;
                isjMin = j == 0;
                
                i1 = isiMax ? 0 : i + 1;
                j1 = isjMax ? 0 : j + 1;
                i2 = isiMin ? iMax : i - 1;
                j2 = isjMin ? jMax : j - 1;
                
                double position[14] = {
                    
                    x[(i * netSize + j) * 2],       // Node (i, j)
                    x[(i * netSize + j) * 2 + 1],
                    x[(i * netSize + j1) * 2],      // Node (i, j+1)
                    x[(i * netSize + j1) * 2 + 1],
                    x[(i1 * netSize + j) * 2],      // Node (i+1, j)
                    x[(i1 * netSize + j) * 2 + 1],
                    x[(i1 * netSize + j2) * 2],     // Node (i+1, j-1)
                    x[(i1 * netSize + j2) * 2 + 1],
                    x[(i * netSize + j2) * 2],      // Node (i, j-1)
                    x[(i * netSize + j2) * 2 + 1],
                    x[(i2 * netSize + j) * 2],      // Node (i-1, j)
                    x[(i2 * netSize + j) * 2 + 1],
                    x[(i2 * netSize + j1) * 2],     // Node (i-1, j+1)
                    x[(i2 * netSize + j1) * 2 + 1],
                
                };
                
                double springs[6] = {
                    
                    spr[i][j][0],
                    spr[i][j][1],
                    spr[i][j][2],
                    spr[i][j2][0],
                    spr[i2][j][1],
                    spr[i2][j1][2]
                
                };
                
                for (int k = 1; k < 7; k++) {
                
                    double prefactor = 0.5 * springs[k - 1] / RESTLEN;
                    current += prefactor * deltaLSqrd(position, i, j, k);
                    position[0] += DEL;
                    dxtemp += prefactor * deltaLSqrd(position, i, j, k);
                    position[0] -= DEL;
                    position[1] += DEL;
                    dytemp += prefactor * deltaLSqrd(position, i, j, k);
                    position[1] -= DEL;
                }
                
                dx[(i * netSize + j) * 2] = (dxtemp - current) / DEL;
                dx[(i * netSize + j) * 2 + 1] = (dytemp - current) / DEL;

            }   
        }
        
    } // End of df();
    
    private:
    
    // deltaLSqrd returns the square change in the length of the spring.
    
    double deltaLSqrd(double * pos, const int row, const int col, const int k) {
        
        double lijk = euclDist(pos, k);
        double l0;
        
        switch (k) {
            case 1: case 2: case 3:
                l0 = rlen[row][col][k - 1];
                break;
            case 4: 
                l0 = rlen[row][j2][0];
                break;
            case 5:
                l0 = rlen[i2][col][1];
                break;
            case 6:
                l0 = rlen[i2][j1][2];
                break;
        }
        
        return (lijk - l0) * (lijk - l0);
        
    }
    
    // euclDistSqrd is a Euclidean distance calculator that works with the
    // pointer nPtr declared in operator(). It takes into account the periodic
    // boundary condition for the system. This is the function that incorporates
    // strain into the network.
    // 
    // Note that we use jMax + 1 and iMax + 1 in these expressions because they
    // are equivalent to the size of the network netSize. 
    
    double euclDist(const double * pos, const int k) {
        
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
                
                    xshift += (jMax + 1) * RESTLEN / 2.0 * (1 + strain * sqrt(3.0));
                    yshift += (iMax + 1) * RESTLEN * sqrt(3.0) / 2.0;
                    
                }
                break;
                
            // Case where nPtr[3] is in row 0 and/or column jMax.
            case 3:
                if (isiMax) {
                
                    xshift += (jMax + 1) * RESTLEN / 2.0 * (1 + strain * sqrt(3.0));
                    yshift += (iMax + 1)* RESTLEN * sqrt(3.0) / 2.0;
                    
                }
                if (isjMin)
                    
                    xshift -= (jMax + 1) * RESTLEN;

                break;
                
            // Case where nPtr[4] is in column jMax.
            case 4:
                if (isjMin) {
                    
                    xshift -= (jMax + 1) * RESTLEN;
                    
                }
                break;
                
            // Case where nPtr[5] is in row iMax.
            case 5: 
                if (isiMin) {
                    
                    yshift -= (iMax + 1) * RESTLEN * sqrt(3.0) / 2.0;
                    xshift -= (jMax + 1) * RESTLEN / 2.0 * (1 + strain * sqrt(3.0));
                    
                }
                break;

            // Case where nPtr[6] is in row iMax and/or column 0.
            case 6:
                if (isiMin) {
                
                    yshift -= (iMax + 1) * RESTLEN * sqrt(3.0) / 2.0;
                    xshift -= (jMax + 1) * RESTLEN / 2.0 * (1 + strain * sqrt(3.0));
                    
                }
                if (isjMax)

                    xshift += (jMax + 1) * RESTLEN;
                
                break;
        }
        
        // xadj and yadj refer to the nodes around pos[0/1].
        
        return sqrt((xadj - pos[0]) * (xadj - pos[0]) + (yadj - pos[1]) * (yadj - pos[1]));
    }
};

#endif /*FUNCD_H_*/
