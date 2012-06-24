/* batch.cpp
 * --------
 *
 * Author: Miles Yucht
 * Date: Mon June 23 2012
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
#include "nonaffinity.h"

using namespace std;

int netSize;
double strain;

int main (int argc, char *argv[]) {

    // Default values.

    double pBond = 0.8;
    double youngMod = 1.0;
    netSize = 20;
    strain = 0.0;

    srand( time(NULL) );
    
    for (; pBond > 0.5; pBond -= 0.025) {
        for (strain = 0.01; strain < 0.09; strain += 0.03) {

            // Now that those are parsed, we can start to generate our network.

            double * position = new double [2 * netSize * netSize];
            double *** sprstiff = new double ** [netSize];
            double *** restlen  = new double ** [netSize];

            for (int i = 0; i < netSize; i++) {

                sprstiff[i] = new double * [netSize];
                restlen[i]  = new double * [netSize];

                for (int j = 0; j < netSize; j++) {

                    // x-coordinate
                    position[(i * netSize + j) * 2] = RESTLEN * (i / 2.0 + j);
                    // predicts initial positoin from strain
                    position[(i * netSize + j) * 2] += i * strain * sqrt(3.0) / 2.0;

                    // y-coordinate
                    position[(i * netSize + j) * 2 + 1] = sqrt(3) / 2 * RESTLEN * i;

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
            printf("%f, %f\n", pBond, strain);

            // Print out position vector to "position_data.txt" in the following format:
            // row,column,xval,yval

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

                        for (int k = 0; k < 3; k++)
                            posFile << "," << restlen[i][j][k];

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
        }
    }

    return 0;
}
