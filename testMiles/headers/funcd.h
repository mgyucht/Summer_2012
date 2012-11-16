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
#include "utils.h"

extern const double RESTLEN;

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
    int i1, i2, j1, j2;
    
    // Constructor. Taking the correct approach.

    Funcd(double *** ssprstiff) : spr(ssprstiff) {}

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
                
                for (int k = 1; k < 4; k++) {
                    
                    double delta = deltaL(pos, k);
                    double temp = 0.5 * spr[i][j][k - 1] * delta * delta;
                    funcvalue += temp;
                
                }
                    
            } // Cycling all columns in a row
        } // Cycling all rows
        
        return funcvalue;
    }

    // df(double *x, double *dx) is the gradient function for the
    // network. This determines the gradient function for a point in the
    // network.
    
    void df(double* x, double* dx) {
    
        iMax = netSize - 1;
        jMax = netSize - 1;
        
        for (int i = 0; i <= iMax; i++) {
            for (int j = 0; j <= jMax; j++) {
                
                double gradient_x = 0.0, gradient_y = 0.0;
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
                    
                    double temp_pos_x = position[2 * k] + xshift(k);
                    double temp_pos_y = position[2 * k + 1] + yshift(k);
                    
                    gradient_x += -springs[k - 1] * deltaL(position, k) 
                        / euclDist(position, k) * (temp_pos_x - position[0]);
                    
                    gradient_y += -springs[k - 1] * deltaL(position, k) 
                        / euclDist(position, k) * (temp_pos_y - position[1]);
                }
                
                dx[(i * netSize + j) * 2] = gradient_x;
                dx[(i * netSize + j) * 2 + 1] = gradient_y;
                
            }
        }
        
    }
    
    private:
    
    // deltaL returns the change in the length of the spring.
    
    double deltaL(double * pos, const int k) {
        
        return euclDist(pos, k) - RESTLEN;
        
    }
    
    // euclDistSqrd is a Euclidean distance calculator that works with the
    // pointer nPtr declared in operator(). It takes into account the periodic
    // boundary condition for the system. This is the function that incorporates
    // strain into the network.
    
    double euclDist(const double * pos, const int k) {
        
        double xadj = pos[2 * k] + xshift(k); 
        double yadj = pos[2 * k + 1] + yshift(k);
        
        return sqrt((xadj - pos[0]) * (xadj - pos[0]) + (yadj - pos[1]) 
                * (yadj - pos[1]));
        
    }

    double xshift(const int &k) {

        double xshift = 0.0;

        switch (k) {

            // Case where nPtr[1] is in column 0.
            case 1: 
                if (isjMax)

                    xshift += (jMax + 1) * RESTLEN;

                break;

                // Case where nPtr[2] is in row 0.
            case 2: 
                if (isiMax) {

                    xshift += RESTLEN * (jMax + 1) / 2.0 * (1 + sqrt(3.0) * strain);

                }
                break;

                // Case where nPtr[3] is in row 0 and/or column jMax.
            case 3:
                if (isiMax) {

                    xshift += RESTLEN * (jMax + 1) / 2.0 * (1 + sqrt(3.0) * strain);

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

                    xshift -= RESTLEN * (jMax + 1) / 2.0 * (1 + sqrt(3.0) * strain);                    

                }
                break;

                // Case where nPtr[6] is in row iMax and/or column 0.
            case 6:
                if (isiMin) {

                    xshift -= RESTLEN * (jMax + 1) / 2.0 * (1 + sqrt(3.0) * strain);                    

                }
                if (isjMax)

                    xshift += (jMax + 1) * RESTLEN;

                break;
        }

        return xshift;
    }

    double yshift(const int &k) {

        double yshift = 0.0;

        switch (k) {

            case 2: 
            case 3:
                if (isiMax) {

                    yshift += (iMax + 1) * RESTLEN * sqrt(3.0) / 2.0;

                }
                break;


                // Case where nPtr[5] is in row iMax.
            case 5: 
            case 6:
                if (isiMin) {

                    yshift -= (iMax + 1) * RESTLEN * sqrt(3.0) / 2.0;

                }
                break;

            default: 
                yshift = 0.0;
                break;
        }

        return yshift;
    }

};

#endif /*FUNCD_H_*/
