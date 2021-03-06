// network.cpp
// ---------
//
// network.cpp contains the data structure for integrate.cpp. These include
// the force calculator, the stress calculator, and the node mover.

#include "network.h"
#include <iostream>
#include <math.h>

// euclDist is a Euclidean distance calculator that works with the
// pointer nPtr declared in operator(). It takes into account the periodic
// boundary condition for the system. This is the function that incorporates
// strain into the network.
// 
// Note that we use jMax + 1 and iMax + 1 in these expressions because they
// are equivalent to the size of the network netSize. 

double euclDist(const double* pos, const int k) {

    double x = (pos[2 * k] - pos[0]);
    double y = (pos[2 * k + 1] - pos[1]);

    return sqrt(x * x + y * y);

}

// deltaLSqrd returns the square change in the length of the spring.

double deltaL(double* pos, const int k) {

    return euclDist(pos, k) - RESTLEN;

}

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

                double delta = deltaL(tempPos, k);
                funcvalue = 0.5 * spring[i][j][k - 1] / RESTLEN * delta * delta;

            }

        } // Cycling all columns in a row
    } // Cycling all rows

    return funcvalue;

}

void Network::getNetForces(Motors motorarray) {

    motorarray.step_motors();

    for (int i = 0; i <= iMax; i++) {

        for (int j = 0; j <= jMax; j++) {

            isiMin = i == 0;
            isiMax = i == iMax;
            isjMin = j == 0;
            isjMax = j == jMax;

            int j1 = isjMax ? 0 : j + 1;
            int i1 = isiMax ? 0 : i + 1;
            int j2 = isjMin ? jMax : j - 1;

            double tempPos[8] = {

                pos[(i * netSize + j) * 2], 
                pos[(i * netSize + j) * 2 + 1], 
                pos[(i * netSize + j1) * 2],
                pos[(i * netSize + j1) * 2 + 1] ,
                pos[(i1 * netSize + j) * 2],
                pos[(i1 * netSize + j) * 2 + 1],
                pos[(i1 * netSize + j2) * 2],
                pos[(i1 * netSize + j2) * 2 + 1]

            };

            // Calculate the net x and y force on each node. Similar to gradient function.

            for (int k = 1; k < 4; k++) {

                double x_displacement = tempPos[2 * k] + xshift(k) - tempPos[0];
                double y_displacement = tempPos[2 * k + 1] + yshift(k) - tempPos[1];

                double dist = euclDist(tempPos, k);

                double cosx = x_displacement / dist;
                double sinx = y_displacement / dist;

                double motorforce = motorarray.getforce(i, j, k);
                double temp = spring[i][j][k - 1] * deltaL(tempPos, k)
                    / RESTLEN + motorforce;

                double xcomp = temp * cosx;
                double ycomp = temp * sinx;
                forces[i][j][2 * k - 2] = xcomp > 1e-10 ? xcomp : 0;
                forces[i][j][2 * k - 1] = ycomp > 1e-10 ? ycomp : 0;

            }

        }

    }

}

void Network::getNetForces() {

    for (int i = 0; i <= iMax; i++) {

        for (int j = 0; j <= jMax; j++) {

            isiMin = i == 0;
            isiMax = i == iMax;
            isjMin = j == 0;
            isjMax = j == jMax;

            int j1 = isjMax ? 0 : j + 1;
            int i1 = isiMax ? 0 : i + 1;
            int j2 = isjMin ? jMax : j - 1;

            double tempPos[8] = {

                pos[(i * netSize + j) * 2], 
                pos[(i * netSize + j) * 2 + 1], 
                pos[(i * netSize + j1) * 2] + xshift(1),
                pos[(i * netSize + j1) * 2 + 1] + yshift(1),
                pos[(i1 * netSize + j) * 2] + xshift(2),
                pos[(i1 * netSize + j) * 2 + 1] + yshift(2),
                pos[(i1 * netSize + j2) * 2] + xshift(3),
                pos[(i1 * netSize + j2) * 2 + 1] + yshift(3)

            };

            // Calculate the net x and y force on each node. Similar to gradient function.

            for (int k = 1; k < 4; k++) {

                double x_displacement = tempPos[2 * k] - tempPos[0];
                double y_displacement = tempPos[2 * k + 1] - tempPos[1];

                double dist = euclDist(tempPos, k);

                double cosx = x_displacement / dist;
                double sinx = y_displacement / dist;

                double temp = spring[i][j][k - 1] * deltaL(tempPos, k)
                    / RESTLEN;

                forces[i][j][2 * k - 2] = temp * cosx;
                forces[i][j][2 * k - 1] = temp * sinx;
            }

        }

    }

}

double Network::calcStress(double strain_rate) {

    double stress = 0.0;
    double prefactor = 1 / (sqrt(3.0) / 2.0 * netSize * netSize);
    double xforce, ydist;

    for (int i = 0; i <= iMax; i++) {

        for (int j = 0; j <= jMax; j++) {

            isiMax = i == iMax;
            isjMin = j == 0;
            isjMax = j == jMax;

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
                ydist = tempPos[2 * k + 1] + yshift(k) - tempPos[1];

                stress += xforce * ydist;

            }

        }

    }

    // Add in viscous network deformation strain.

    stress = stress * prefactor; //  + ETA * strain_rate;

    return stress;

}

double affvx(int r, double s_rate)
{
    double hmid = (netSize-1.0)/2.0;
    //return sqrt(3.0) / 2.0 * s_rate * ((r-hmid) / hmid);
    return sqrt(3.0) / 2.0 * s_rate * (r-hmid);
}

void Network::moveNodes(double shear_rate, double temp) {

    double netx, nety, affvel, gamma;
    double d = KB * temp / (6 * PI * ETA * RADIUS);
    double sigma = sqrt(2 * d * TIMESTEP);

    affvel = affvx(netSize - 1, shear_rate);
    affdel += affvel * TIMESTEP;

    for (int i = 0; i <= iMax; i++) {

        for (int j = 0; j <= jMax; j++) {

            int currentx = (i * netSize + j) * 2;
            int currenty = currentx + 1;

            isjMax = j == jMax;
            isiMin = i == 0;
            isjMin = j == 0;

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

            netx = fHooke[0] + fHooke[2] + fHooke[4] - fHooke[6] - fHooke[8] - fHooke[10];
            nety = fHooke[1] + fHooke[3] + fHooke[5] - fHooke[7] - fHooke[9] - fHooke[11];

            affvel = affvx(i, shear_rate);
            // vel_fluid = sqrt(3.0) / 4.0 * netSize * shear_rate * (2 * ((double) i - netSize) / (netSize + 1) + 1);

            gamma = 4 * PI * ETA * RADIUS;

            if (temp > 1e-15)
            {
                double theta = 2 * PI * randDouble(0, 1);
                double r = sigma * sqrt(-2 * log(randDouble(0, 1)));

                double temp_fluct_x = r * cos(theta);
                double temp_fluct_y = r * sin(theta);

                delta[currentx] = TIMESTEP * (netx / gamma + affvel) + temp_fluct_x;
                delta[currenty] = TIMESTEP * (nety / gamma) + temp_fluct_y;
                pos[currentx] += delta[currentx];
                pos[currenty] += delta[currenty];
            } else // temp = 0.0
            {
                delta[currentx] = TIMESTEP * (netx / gamma + affvel);
                delta[currenty] = TIMESTEP * (nety / gamma);
                pos[currentx] += delta[currentx];
                pos[currenty] += delta[currenty];
            }

            if (pos[currentx] != pos[currentx] || pos[currenty] != pos[currenty])
            {
                throw("NaN value assigned");
            }

        }

    }

}

double Network::xshift(const int &k) {

    double xshift = 0.0;
    double netwidth = netSize;
    double d = affdel; // delta for the entire network
    double netshift = netSize / 2.0 + (2.0 + 2.0 / (netSize - 1.0)) * d;

    switch (k) {

        // Case where nPtr[1] is in column 0.
        case 1: 
            if (isjMax)

                xshift += netwidth;

            break;

        // Case where nPtr[2] is in row 0.
        case 2: 
            if (isiMax) {

                xshift += netshift;

            }
            break;

        // Case where nPtr[3] is in row 0 and/or column jMax.
        case 3:
            if (isiMax) {

                xshift += netshift;

            }
            if (isjMin)

                xshift -= netwidth;

            break;

        // Case where nPtr[4] is in column jMax.
        case 4:
            if (isjMin) {

                xshift -= netSize;

            }
            break;

        // Case where nPtr[5] is in row iMax.
        case 5: 
            if (isiMin) {

                xshift -= netshift;

            }
            break;

        // Case where nPtr[6] is in row iMax and/or column 0.
        case 6:
            if (isiMin) {

                xshift -= netshift;

            }
            if (isjMax)

                xshift += netwidth;

            break;
    }

    return xshift;
}

double Network::yshift(const int &k) {

    double yshift = 0.0;

    switch (k) {

        case 2: 
        case 3:
            if (isiMax) {

                yshift += netSize * sqrt(3.0) / 2.0;

            }
            break;

        // Case where nPtr[5] is in row iMax.
        case 5: 
        case 6:
            if (isiMin) {

                yshift -= netSize * sqrt(3.0) / 2.0;

            }
            break;

        default: 
            break;
    }

    return yshift;
}
