// nonaffinity.cpp
// ---------------
//
// nonaffinity.cpp takes as its input a position file and outputs the nonaffinity measure 
// for the network. The nonaffinity measure of a network is a way of measuring how much 
// a network deforming nonaffinely varies from the affine deformation. It is calculated 
// by the formula:
//
//          1
//     Γ = --- Σ (u_i - uaff_i)²
//         Nγ²
//  
// where Γ is the nonaffinity measure, N is the number of nodes in the network, γ is the
// strain on the network, u_i is the displacement of the node in the nonaffine 
// deformation, and uaff_i the displacement ofthe same node in the affine one.

#include <cstdlib>
#include <fstream>
#include <sstream>
#include <cmath>
#include <cstring>

#include "nonaffinity.h"

double nonAffinity(double* position, int netSize, double strain)
{
    int row, col;
    double xval, yval, currentx, currenty, prefactor, sqrdisp = 0, 
           nonaffinity = 0;

    if (std::abs(strain) < 1E-15) 
    {
        return 0.0;
    }
    
    prefactor = 1 / (netSize * netSize * strain * strain);
    
    for (row = 0; row < netSize; row++)
    {
        for (col = 0; col < netSize; col++)
        {
            currentx = position[(row * netSize + col) * 2];
            currenty = position[(row * netSize + col) * 2 + 1];
            
            xval = (row / 2.0 + col) + row * sqrt(3.0) * strain / 2.0;
            yval = sqrt(3.0) / 2.0 * row;
            
            sqrdisp += (currentx - xval) * (currentx - xval) + (currenty - yval)
                            * (currenty - yval);
        }
    }
    
    nonaffinity = prefactor * sqrdisp;
    
    return nonaffinity;
    
}
