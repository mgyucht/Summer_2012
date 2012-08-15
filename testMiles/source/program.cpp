/* program.cpp
 * --------
 *
 * program.cpp is the client file in simulating the spring networks.
 * Usage: program -str <strain> -size <network size> -p <bond probability> -y \
 * <young's modulus for springs>
 *
 * Author: Miles Yucht
 * Date: Mon June 27 2012
 */

#include <string>
#include <ctime>
#include <cstdio>
#include <time.h>
#include <stdlib.h>
#include <math.h>
#include "funcd.h"
#include "frprmn.h"
#include "utils.h"
#include "network.h"
#include "print.h"

using namespace std;

const double RESTLEN = 1;

int netSize;
double strain;

int main (int argc, char *argv[]) {

    // Default values.

    double pBond = 0.8;
    double youngMod = 1.0;
    netSize = 20;
    strain = 0.0;

    srand( time(NULL) );

    //Otherwise, parse the command line parameters.

    for (int i = 1; i < argc; i += 2) {
        string str = argv[i];
        if (!str.compare("-str")) {
            strain = atof(argv[i + 1]);
        } else if (!str.compare("-size")) {
            netSize = atoi(argv[i + 1]);
        } else if (!str.compare("-p")) {
            pBond = atof(argv[i + 1]);
        } else if (!str.compare("-y")) {
            youngMod = atof(argv[i + 1]);
        } else if (!str.compare("--help") || !str.compare("help")) {
            usageExit();
        } else {
            printf("%s is an illegal argument.\n", argv[i]);
            usageExit();
        }
    }

    printf("strain is %3.2f, pBond is %4.3f\n", strain, pBond);
    // Now that those are parsed, we can start to generate our network.

    double network_height = netSize * sqrt(3.0) / 2.0;
    double strain_distance = strain * network_height;
    
    double position[2 * netSize * netSize];
    double *** sprstiff = new double ** [netSize];
    double *** forces = new double ** [netSize];

    for (int i = 0; i < netSize; i++) {

        sprstiff[i] = new double * [netSize];
        forces[i] = new double * [netSize];

        for (int j = 0; j < netSize; j++) {

            // x-coordinate
            position[(i * netSize + j) * 2] = RESTLEN * (i / 2.0 + j);
            // predicts initial position from strain
            position[(i * netSize + j) * 2] += strain_distance * i / (double) netSize;

            // y-coordinate
            position[(i * netSize + j) * 2 + 1] = sqrt(3) / 2 * RESTLEN * i;

            sprstiff[i][j] = stiffVecGen(pBond, youngMod, 3);
            forces[i][j] = new double[6];
        }
    }
    
    // Minimize energy

    Funcd test(sprstiff);
    Frprmn<Funcd> frprmn(test);
    double *newArray = frprmn.minimize(position);
    double newEnergy = test(newArray);
    
    Network myNetwork(position, sprstiff, sprstiff, forces, 0.01);
    myNetwork.getNetForces();
    double stress = myNetwork.calcStress();
    
    FILE * compareData = fopen("compare_data.txt", "a");
    
    fprintf(compareData, "Stress                 = %f\n", stress);
    fprintf(compareData, "G (stress calculation) = %f\n", stress / strain);
    fprintf(compareData, "G (energy calculation) = %f\n", 2 * newEnergy / 
            (sqrt(3.0) / 2.0 * netSize * netSize * strain * strain));
    fflush(compareData);
    
    fclose(compareData);
    

    // Print out position vector to "position_data.txt" in the following format:
    // row,column,xval,yval

    Printer myPrinter(myNetwork, youngMod, pBond);
    myPrinter.printPos("position_data.txt");
    myPrinter.printNonAff("nonaff_data.txt");
    myPrinter.printEnergy("energy_data.txt", newEnergy);
    myPrinter.printStress("stress_data.txt", &stress, 1);

    // Cleanup

    for (int i = 0; i < netSize; i++) {
        for (int j = 0; j < netSize; j++) {
            delete[] sprstiff[i][j];
        }
        delete[] sprstiff[i];
    }

    delete[] sprstiff;

    return 0;
}
