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

#include <iomanip>
#include <string>
#include <ctime>
#include <time.h>
#include <stdlib.h>
#include <math.h>
#include <sstream>
#include "compileinfo.h"
#include "utils.h"
#include "network.h"
#include "nonaffinity.h"
#include "print.h"

using namespace std;

extern const string output_path;

// Rest length for springs.
extern const double RESTLEN;
// Viscosity of the fluid.
extern const double ETA;
// Radius for Stokes' drag.
extern const double RADIUS;
// Young's modulus for springs.
extern const double YOUNGMOD;
// Time step for the simulation.
double TIMESTEP;

// NOTE: One unit of time in this simulation is equivalent to 10^-5 seconds in
// reality: 1 s^* = 10^-5 s.

int netSize;
double strain;

int main (int argc, char *argv[]) {

    // Default values.

    int prngseed = 0, nTimeSteps = 20000, job;
    double pBond = 0.8, strRate = 1.0, currentTime = 0.0, initStrain = 0.01;
    
    netSize = 20;
    
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
            
        } else if (!str.compare("-j")) {
            
            job = atoi(argv[i + 1]);
            
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
    
    TIMESTEP = 1 / (1000 * strRate);
    
    // Get filenames ready.
    
    ostringstream convert;
    convert << job << "/";
    
    string root_path = output_path + convert.str();
    
    // Make the directory if it doesn't exist.
    
    string posFilePath    = root_path + posFileName + extension;
    string nonaffFilePath = root_path + nonaffFileName + extension;
    string stressFilePath = root_path + stressFileName + extension;
    string energyFilePath = root_path + energyFileName + extension;
    
    char* posFileFull    = (char *) (posFilePath).c_str();
    char* nonaffFileFull = (char *) (nonaffFilePath).c_str();
    char* stressFileFull = (char *) (stressFilePath).c_str();
    char* energyFileFull = (char *) (energyFilePath).c_str();
    
    // Initialize PRNG.
    if (prngseed == 0) 

        srand( time(NULL) );
    
    else 

        srand( prngseed );

    // Now that those are parsed, we can start to generate our network.
     
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

            // y-coordinate
            position[(i * netSize + j) * 2 + 1] = sqrt(3) / 2 * RESTLEN * i;

            sprstiff[i][j] = stiffVecGen(pBond, 3);
            velocities[i][j] = new double[2];
            netForces[i][j] = new double[6];
            
        }
    }
    
    double* strain_array = new double[nTimeSteps];
    double* strain_rate = new double[nTimeSteps];
    
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
    
    // myPrinter.printPos(posFileFull);
    // myPrinter.printNonAff(nonaffFileFull);
    myPrinter.printStress(stressFileFull, stress_array, strain_array);
    // myPrinter.printEnergy(energyFileFull, newEnergy);

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
    delete[] strain_rate;
    delete[] strain_array;

    return 0;
}
