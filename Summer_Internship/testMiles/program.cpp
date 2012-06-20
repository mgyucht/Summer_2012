/* program.cpp
 * --------
 *
 * program.cpp is the client file in simulating the spring networks.
 * Usage: program -str <strain> -size <network size> -p <bond probability> -y \
 * <young's modulus for springs>
 *
 * Author: Miles Yucht
 * Date: Mon June 11 2012
 */

/*
 * Working revisions:
 *
 * Only two dimensions.
 * nr3.h is no longer needed! Yay!
 *
 */

#include <string>
#include <ctime>
#include <time.h>
#include <stdlib.h>
#include <math.h>
#include <fstream>
#include "funcd.h"
#include "frprmn.h"
#include "utils.h"
#include "debug.h"

using namespace std;

int netSize;
double strain;

int main (int argc, char *argv[]) {

    // Default values.

    double pBond = 1.0;
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
        } else if (!str.compare("-help") || !str.compare("help")) {
            usageExit();
        } else {
            printf("%s is an illegal argument.\n", argv[i]);
            usageExit();
        }
    }

    // Now that those are parsed, we can start to generate our network.

    double * position = new double [2 * netSize * netSize];
    double *** sprstiff = new double ** [netSize];
    double *** restlen  = new double ** [netSize];

    for (int i = 0; i < netSize; i++) {

        sprstiff[i] = new double * [netSize];
        restlen[i]  = new double * [netSize];

        for (int j = 0; j < netSize; j++) {

            // x-coordinate
            position[i * 2 * netSize + j * 2] = RESTLEN * (i / 2.0 + j);

            // y-coordinate
            position[i * 2 * netSize + j * 2 + 1] = sqrt(3) / 2 * RESTLEN * i;

            sprstiff[i][j] = stiffVecGen(pBond, youngMod, 3);

            restlen[i][j] = new double [3];

            for (int k = 0; k < 3; k++)
                restlen[i][j][k] = RESTLEN;
        }
    }

    // Minimize energy

    Funcd test(sprstiff, restlen);
    Frprmn<Funcd> frprmn(test);
    double *newArray = frprmn.minimize(position);
    double newEnergy = test(newArray);
    printf("E = %f\n", newEnergy);

    // Print out position vector to "position_data.txt" in the following format:
    // row,column,xval,yval

#if BATCH

    ofstream posFile("position_data.txt", ios::app);

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

                for (int k = 0; k < 3; k++)
                    posFile << "," << restlen[i][j][k];

                posFile << endl;
            }

        }

    }

    posFile.close();

    // Print out energy and p for the experiment

    ofstream engFile("energy_data.txt", ios::app);

    if (engFile.is_open()) {

        engFile << newEnergy << "," << pBond << "," << strain << endl;

    }

    engFile.close();

#endif

    // Cleanup

    delete[] position;

    for (int i = 0; i < netSize; i++) {
        for (int j = 0; j < netSize; j++) {
            delete[] sprstiff[i][j];
            delete[] restlen[i][j];
        }
        delete[] sprstiff[i];
        delete[] restlen[i];
    }

    delete[] sprstiff;
    delete[] restlen;

    return 0;
}
