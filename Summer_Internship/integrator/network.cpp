// funcs.cpp
// ---------
//
// funcs.cpp contains the essential functions for integrate.cpp. These include
// the force calculator, the stress calculator, and the node mover.

#include <math.h>
#include "network.h" 

// Network Function Definitions.
double Network::operator() () {
    
    double funcvalue = 0;

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

            int j1 = isjMax ? 0 : j + 1;
            int i1 = isiMax ? 0 : i + 1;
            int j2 = isjMin ? jMax : j - 1;

            double tempPos[8] = {

                pos[(i * netSize + j) * 2], 
                pos[(i * netSize + j) * 2 + 1], 
                pos[(i * netSize + j1) * 2],
                pos[(i * netSize + j1) * 2 + 1],
                pos[(i1 * netSize + j) * 2],
                pos[(i1 * netSize + j) * 2 + 1],
                pos[(i1 * netSize + j2) * 2],
                pos[(i1 * netSize + j2) * 2 + 1]

            };

            for (int k = 1; k < 4; k++) {

                funcvalue = 0.5 * spr[i][j][k - 1] / RESTLEN * deltaLSqrd(tempPos, i, j, k);

            }

        } // Cycling all columns in a row
    } // Cycling all rows

    return funcvalue;

}

void Network::getNetForces() {

    for (int i = 0; i < iMax; i++) {

        for (int j = 0; j < jMax; j++) {
            
            isiMin = i == 0;
            isiMax = i == iMax;
            isjMin = j == 0;

            int j1 = isjMax ? 0 : j + 1;
            int i1 = isiMax ? 0 : i + 1;
            int j2 = isjMin ? jMax : j - 1;

            double tempPos[8] = {

                pos[(i * netSize + j) * 2], 
                pos[(i * netSize + j) * 2 + 1], 
                pos[(i * netSize + j1) * 2],
                pos[(i * netSize + j1) * 2 + 1],
                pos[(i1 * netSize + j) * 2],
                pos[(i1 * netSize + j) * 2 + 1],
                pos[(i1 * netSize + j2) * 2],
                pos[(i1 * netSize + j2) * 2 + 1]

            };

            // Calculate the net x and y force on each node. Similar to gradient function.

            for (int k = 1; k < 4; k++) {

                forces[i][j][2 * k - 2] = spr[i][j][k - 1] 
                    * cos( atan2( tempPos[2 * k + 1] - tempPos[1], tempPos[2 * k] - tempPos[0] ) )
                    * sqrt(deltaLSqrd(tempPos, i, j, k));

                forces[i][j][2 * k - 1] = spr[i][j][k - 1] 
                    * sin( atan2( tempPos[2 * k + 1] - tempPos[1], tempPos[2 * k] - tempPos[0] ) )
                    * sqrt(deltaLSqrd(tempPos, i, j, k));

            }

        }

    }

}

double Network::calcStress() {

    double stress = 0.0;
    double prefactor = 1 / (sqrt(3.0) / 2.0 * netSize * netSize);
    double xforce, ydist;

    for (int i = 0; i < iMax; i++) {

        for (int j = 0; j < jMax; j++) {

            isiMin = i == 0;
            isiMax = i == iMax;
            isjMin = j == 0;
            
            int j1 = isjMax ? 0 : j + 1;
            int i1 = isiMax ? 0 : i + 1;
            int j2 = isjMin ? jMax : j - 1;

            double tempPos[8] = {

                pos[(i * netSize + j) * 2], 
                pos[(i * netSize + j) * 2 + 1], 
                pos[(i * netSize + j1) * 2],
                pos[(i * netSize + j1) * 2 + 1],
                pos[(i1 * netSize + j) * 2],
                pos[(i1 * netSize + j) * 2 + 1],
                pos[(i1 * netSize + j2) * 2],
                pos[(i1 * netSize + j2) * 2 + 1]

            };

            for (int k = 1; k < 4; k++) {

                // Get the x-component of the force.
                xforce = forces[i][j][2 * k - 2];

                // Get the y-distance between nodes.
                ydist = tempPos[2 * k + 1] - tempPos[1];

                stress += xforce * ydist;

            }

        }

    }

    return stress * prefactor;

}

void Network::moveNodes() {

    double netx, nety;

    for (int i = 0; i < iMax; i++) {

        for (int j = 0; j < jMax; j++) {

            isiMin = i == 0;
            isiMax = i == iMax;
            isjMin = j == 0;
            isjMax = j == jMax;
            
            int i1 = isiMax ? 0 : i + 1;
            int j1 = isjMax ? 0 : j + 1;
            int i2 = isiMin ? iMax : i - 1;
            int j2 = isjMin ? jMax : j - 1;

            double fHooke[12] = {

                forces[i][j][0],
                forces[i][j][1],
                forces[i][j][2],
                forces[i][j][3],
                forces[i][j][4],
                forces[i][j][5],
                forces[i][j2][0],
                forces[i][j2][1],
                forces[i2][j][2],
                forces[i2][j][3],
                forces[i2][j1][4],
                forces[i2][j1][5]

            };

            netx = fHooke[0] + fHooke[2] + fHooke[4] + fHooke[6] + fHooke[8] + fHooke[10];
            nety = fHooke[1] + fHooke[3] + fHooke[5] + fHooke[7] + fHooke[9] + fHooke[11];

            pos[(i * netSize + j) * 2] += vels[i][j][0] * timestep;
            pos[(i * netSize + j) * 2 + 1] += vels[i][j][1] * timestep;

            vels[i][j][0] += netx / MASS * timestep;
            vels[i][j][1] += nety / MASS * timestep;

        }

    }

}

// euclDistSqrd is a Euclidean distance calculator that works with the
// pointer nPtr declared in operator(). It takes into account the periodic
// boundary condition for the system. This is the function that incorporates
// strain into the network.
// 
// Note that we use jMax + 1 and iMax + 1 in these expressions because they
// are equivalent to the size of the network netSize. 

double Network::euclDist(const double* pos, const int k) {

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
