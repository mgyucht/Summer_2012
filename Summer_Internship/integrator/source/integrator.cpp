/* integrate.cpp
 * --------
 *
 * integrator.cpp is the client file in simulating the spring networks. This 
 * program works by integrating the equations of motion over all of the nodes
 * in the network. It also takes into account the viscous force from the 
 * shearing of a liquid.
 *
 * Author: Miles Yucht
 * Date: Mon June 29 2012
 */

#include <string>
#include <ctime>
#include <time.h>
#include <stdlib.h>
#include <math.h>
#include <sstream>
#include "utils.h"
#include "network.h"
#include "nonaffinity.h"
#include "print.h"

using namespace std;

// Rest length for springs.
const double RESTLEN = 1.0;
// Viscosity of the fluid.
const double ETA = 1;
// Radius for Stokes' drag.
const double RADIUS = 0.1;
// Young's modulus for springs.
const double YOUNGMOD = 1.0;
// Time step for the simulation.
const double TIMESTEP = 1e-4;

// NOTE: One unit of time in this simulation is equivalent to 10^-5 seconds in
// reality: 1 s^* = 10^-5 s.

int netSize;
double strain;

int main (int argc, char *argv[]) {

    // Default values.

    int prngseed = 0, nTimeSteps = 100;
    double pBond = 0.8, strRate = 1.0, currentTime = 0.0, initStrain = 0.01;
    
    netSize = 20;
    
    string output_path = "output/";
    string energyFileName = "energy_data";
    string posFileName = "position_data";
    string nonaffFileName = "nonaff_data";
    string stressFileName = "stress_data";
    string extension = ".txt";
     
    // If there are any, parse the command line parameters.

    for (int i = 1; i < argc; i += 2) {
        string str = argv[i];
        
        if (!str.compare("-str")) {
            
            initStrain = atof(argv[i + 1]);
            
        } else if (!str.compare("-size")) {
            
            netSize = atoi(argv[i + 1]);
            
        } else if (!str.compare("-p")) {
            
            pBond = atof(argv[i + 1]);
            
        } else if (!str.compare("-rate")) {
            
            strRate = atof(argv[i + 1]);
            
        } else if (!str.compare("-seed")) {
            
            prngseed = atoi(argv[i + 1]);
            
        } else if (!str.compare("-n")) {

            nTimeSteps = atoi(argv[i + 1]);
            
        } else if (!str.compare("-energy-fn")) {
            
            energyFileName = argv[i + 1];
            
        } else if (!str.compare("-position-fn")) {
            
            posFileName = argv[i + 1];
            
        } else if (!str.compare("-nonaff-fn")) {
            
            nonaffFileName = argv[i + 1];
            
        } else if (!str.compare("-stress-fn")) {
            
            stressFileName = argv[i + 1];
            
        } else if (!str.compare("-help") || !str.compare("help")) {
            
            usageExit();
            
        } else {
            
            printf("%s is an illegal argument.\n", argv[i]);
            usageExit();
            
        }
    }
    
    // Initialize PRNG.
    if (prngseed == 0) 

        srand( time(NULL) );
    
    else 

        srand( prngseed );

    // Now that those are parsed, we can start to generate our network.
     
    double network_height = sqrt(3.0) / 2.0 * netSize * netSize;
    double strain_dist = initStrain * network_height;

    double* position = new double [2 * netSize * netSize];
    double* stress_array = new double [nTimeSteps];
    double*** sprstiff = new double** [netSize];
    double*** velocities = new double** [netSize];
    double*** netForces = new double** [netSize];

    for (int i = 0; i < netSize; i++) {

        sprstiff[i] = new double* [netSize];
        velocities[i] = new double* [netSize];
        netForces[i] = new double* [netSize];

        for (int j = 0; j < netSize; j++) {

            // x-coordinate
            position[(i * netSize + j) * 2] = RESTLEN * (i / 2.0 + j);
            // predicts initial position from strain
            position[(i * netSize + j) * 2] += i * strain_dist / network_height 
                * RESTLEN;

            // y-coordinate
            position[(i * netSize + j) * 2 + 1] = sqrt(3) / 2 * RESTLEN * i;

            sprstiff[i][j] = stiffVecGen(pBond, 3);
            velocities[i][j] = new double[2];
            netForces[i][j] = new double[6];
            
        }
    }
    
    double strain_array[nTimeSteps];
    double strain_rate[nTimeSteps];
    
    for (int i = 0; i < nTimeSteps; i++) {
    
        strain_array[i] = initStrain * sin(strRate * i * TIMESTEP);
        strain_rate[i]   = initStrain * strRate * TIMESTEP * cos(strRate * i * TIMESTEP);
    
    }

    // Integrate motion over the nodes.
    // 
    // Technical note: Stress is a rank 2 tensor. However, because our interest
    // is in the shear modulus of the network, we discard the isotropic
    // elements of the stress. Therefore, because this is a 2-D network, there 
    // is only one real term of interest, because σ_xy = σ_yx. This is the 
    // number that is output by calcStress.
    
    Network myNetwork(position, sprstiff, velocities, netForces);
    Printer myPrinter(myNetwork, pBond, nTimeSteps);
    
    for (int i = 0; i < nTimeSteps; i++) {
    
        currentTime += TIMESTEP;
        strain = strain_array[i];
        
        myNetwork.getNetForces();
        stress_array[i] = myNetwork.calcStress();
        myNetwork.moveNodes(strain_rate[i]);
        
    }

    double newEnergy = myNetwork();
    
    // Print out position vector to "position_data.txt" in the following format:
    // row,column,xval,yval, sprstiffs, restlens.

    char* posFileFull = (char *) (output_path + posFileName + extension).c_str();
    myPrinter.printPos(posFileFull);
    
    char* nonaffFileFull = (char *) (output_path + nonaffFileName + extension).c_str();
    myPrinter.printNonAff(nonaffFileFull);
    
    char* stressFileFull = (char *) (output_path + stressFileName + extension).c_str();
    myPrinter.printStress(stressFileFull, stress_array, strain_array);
    
    char* energyFileFull = (char *) (output_path + energyFileName + extension).c_str();
    myPrinter.printEnergy(energyFileFull, newEnergy);

    // Cleanup

    delete[] position;
    delete[] stress_array;

    for (int i = 0; i < netSize; i++) {
        for (int j = 0; j < netSize; j++) {
            delete[] sprstiff[i][j];
            delete[] velocities[i][j];
            delete[] netForces[i][j];
        }
        delete[] sprstiff[i];
        delete[] velocities[i];
        delete[] netForces[i];
    }

    delete[] sprstiff;
    delete[] velocities;
    delete[] netForces;

    return 0;
}
