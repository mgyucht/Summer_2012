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

#include <stdlib.h>
#include <fstream>
#include <sstream>
#include <math.h>
#include <cstring>

using namespace std;

double nonAffinity(const char * argv) {
    
    int netSize, row, col, count = 0;
    double strain, rlen, xval, yval, currentx, currenty, prefactor, sqrdisp = 0;
    double nonaffinity = 0;
    
    ifstream posFile(argv, ifstream::in);
    
    if (posFile.is_open()) {
        
        string line;
        string token;
        stringstream iss;
        
        // Discard the first line of the file.
        getline(posFile, line);
        getline(posFile, line);
    
        iss << line;
        getline(iss, token, ',');
        netSize = atoi(token.c_str());
        getline(iss, token, ',');
        strain = atof(token.c_str());
        
        if (strain < 1E-15) {
            
            return 0.0;
            
        }
        
        prefactor = 1 / (netSize * netSize * strain * strain);
        
        while (count < netSize * netSize) {
            
            // Clear iss and put the next line in it.
            iss.str("");
            getline(posFile, line);
            iss << line;
            
            getline(iss, token, ',');
            istringstream(token) >> row;
            getline(iss, token, ',');
            istringstream(token) >> col;
            getline(iss, token, ',');
            istringstream(token) >> currentx;
            getline(iss, token, ',');
            istringstream(token) >> currenty;
            getline(iss, token, ',');
            getline(iss, token, ',');
            getline(iss, token, ',');
            getline(iss, token, ',');
            istringstream(token) >> rlen;
            
            xval = rlen * (row / 2.0 + col) + row * sqrt(3.0) * strain / 2.0;
            yval = sqrt(3.0) / 2.0 * rlen * row;
            
            sqrdisp += (currentx - xval) * (currentx - xval) + (currenty - yval)
                            * (currenty - yval);
            ++count;
        
        }
        
        nonaffinity = prefactor * sqrdisp;
        
        posFile.close();
            
    } else {
        
        printf("Something bad has happened, and your position file could not be read.\n");
        printf("Are you sure you entered the filename correctly?");
        
    }
    
    return nonaffinity;
    
}
