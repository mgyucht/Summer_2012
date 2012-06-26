/* integrate.cpp
 * --------
 *
 * integrate.cpp is the client file in simulating the spring networks. This 
 * program works by integrating the equations of motion over all of the nodes
 * in the network. It also takes into account the viscous force from the 
 * shearing of a liquid.
 * 
 * Usage: integrate -str <strain> -size <network size> -p <bond probability> -y \
 * <young's modulus for springs>
 *
 * Author: Miles Yucht
 * Date: Mon June 11 2012
 */

#include <string>
#include <ctime>
#include <time.h>
#include <stdlib.h>
#include <math.h>
#include <fstream>
#include "utils.h"
#include "network.h"
#include "debug.h"
#include "nonaffinity.h"

using namespace std;

// Rest length for springs.
const double RESTLEN = 1.0;

int netSize;
double strain;

int main (int argc, char *argv[]) {

    // Default values.

    int prngseed = 0, nTimeSteps = 100;
    double pBond = 0.8, strRate = 1.0, youngMod = 1.0, currentTime = 0.0, 
           timeStep = 0.001, initStrain = 0.01;
    
    netSize = 20;
     
    // If there are any, parse the command line parameters.

    for (int i = 1; i < argc; i += 2) {
        string str = argv[i];
        
        if (!str.compare("-str")) {
            
            initStrain = atof(argv[i + 1]);
            
        } else if (!str.compare("-size")) {
            
            netSize = atoi(argv[i + 1]);
            
        } else if (!str.compare("-p")) {
            
            pBond = atof(argv[i + 1]);
            
        } else if (!str.compare("-y")) {
            
            youngMod = atof(argv[i + 1]);
            
        } else if (!str.compare("-rate")) {
            
            strRate = atof(argv[i + 1]);
            
        } else if (!str.compare("-seed")) {
            
            prngseed = atoi(argv[i + 1]);
            
        } else if (!str.compare("-step")) {

            nTimeSteps = atoi(argv[i + 1]);

        } else if (!str.compare("-help") || !str.compare("help")) {
            
            usageExit();
            
        } else {
            
            printf("%s is an illegal argument.\n", argv[i]);
            usageExit();
            
        }
    }
    
    // Initialize PRNG.
    if (!prngseed) 

        srand( time(NULL) );
    
    else 

        srand( prngseed );

    // Now that those are parsed, we can start to generate our network.

    double* position = new double [2 * netSize * netSize];
    double* stress = new double [nTimeSteps];
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
            position[(i * netSize + j) * 2] += i * strain * sqrt(3.0) / 2.0;

            // y-coordinate
            position[(i * netSize + j) * 2 + 1] = sqrt(3) / 2 * RESTLEN * i;

            sprstiff[i][j] = stiffVecGen(pBond, youngMod, 3);
            velocities[i][j] = new double[2];
            netForces[i][j] = new double[6];
            
        }
    }

    // Integrate motion over the nodes.
    // 
    // Technical note: Stress is a rank 2 tensor. However, because our interest
    // is in the shear modulus of the network, we discard the isotropic
    // elements of the stress. Therefore, because this is a 2-D network, there 
    // is only one real term of interest, because σ_xy = σ_yx. This is the 
    // number that is output by calcStress.
    
    Network myNetwork(position, sprstiff, velocities, netForces, timeStep);
    
    for (int i = 0; i < nTimeSteps; i++) {
    
        strain = initStrain * cos(strRate * currentTime);
        
        myNetwork.getNetForces();
        stress[i] = myNetwork.calcStress();
        myNetwork.moveNodes();
        
        currentTime += timeStep;
        
    }

    printf("E = %f\n", myNetwork());
    
    // Print out position vector to "position_data.txt" in the following format:
    // row,column,xval,yval, sprstiffs, restlens.

    #if PRINTPOS

    string posFileName = "position_data.txt";
    ofstream posFile(posFileName.c_str(), ios::trunc);

    if (posFile.is_open()) {
 
        posFile << "Network Size,Strain,Young's Modulus,Probability" << endl;
        posFile << netSize << "," << strain << "," << youngMod << ","
            << pBond << endl;

        for (int i = 0; i < netSize; i++) {

            for (int j = 0; j < netSize; j++) {

                // Print row, col, position, sprstiff, rlen information to file.

                posFile << i << "," << j << "," << position[(i * netSize + j) * 2]
                    << "," << position[(i * netSize + j) * 2 + 1];

                for (int k = 0; k < 3; k++) 
                    posFile << "," << sprstiff[i][j][k];

                posFile << endl;
            }

        }

    }

    posFile.close();
    
    #if PRINTNONAFFINITY
    
    string nonaffFileName = "nonaff_data.txt";
    ofstream nonaffFile(nonaffFileName.c_str(), ios::app);

    if (nonaffFile.is_open()) {
        
        double nonaff = nonAffinity(posFileName.c_str());
        printf("Γ = %f\n", nonaff);
        nonaffFile << netSize << "," << strain << "," << pBond << "," 
            << nonaff << endl;
    
    }
    
    nonaffFile.close();

    #endif //PRINTNONAFFINITY
    
    #endif //PRINTPOS
    
    #if PRINTENERGY

    // Print out energy and p for the experiment

    ofstream engFile("energy_data.txt", ios::app);

    if (engFile.is_open()) {

        engFile << newEnergy << "," << pBond << "," << strain << endl;

    }

    engFile.close();

    #endif

    // Cleanup

    delete[] position;
    delete[] stress;

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
