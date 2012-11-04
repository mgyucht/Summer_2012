#ifndef NETWORK_H_
#define NETWORK_H_

// network.h
// -------
//
// network.h is the header for network.cpp, which contains the essential methods 
// for integrate.cpp. It defines the struct Network.

#include <math.h>
#include "utils.h"
#include "motors.h"

extern double TIMESTEP;
extern const double RESTLEN;
extern const double ETA;
extern const double RADIUS;

extern int netSize;
extern double strain;

static const double KB = 1;
static const double PI = 3.1415926535;

struct Network {

    double *pos;
    double *delta;
    double ***spring;
    double ***forces;

    int iMax, jMax;

    bool isiMax, isjMax, isiMin, isjMin;

    Network(double *ppos, double *ddelta, double ***sspring, double ***fforces) :
        pos(ppos),
        delta(ddelta),
        spring(sspring),
        forces(fforces),
        isiMax(false), isjMax(false), isiMin(true), isjMin(true) {

        iMax = netSize - 1;
        jMax = netSize - 1;

    }

    double operator() ();

    // getNetForces sets the forces array. For each node at lattice index (i, j),
    // the forces exerted on the node by the nodes (i, j+1), (i+1, j), and 
    // (i+1, j-1) are associated with the node (i, j) in forces by the following
    // table:
    //
    // forces[i][j][0] - x component of force with node (i, j+1)
    // forces[i][j][1] - y component of force with node (i, j+1)
    // forces[i][j][2] - x component of force with node (i+1, j)
    // forces[i][j][3] - y component of force with node (i+1, j)
    // forces[i][j][4] - x component of force with node (i+1, j-1)
    // forces[i][j][5] - y component of force with node (i+1, j-1)
    // 
    // The force is calculated using Hooke's law.

    void getNetForces(Motors /* Motors object */);
    void getNetForces();

    double calcStress(double strain_rate);

    void moveNodes(double shear_rate, double temp);

    private:

    // deltaLSqrd returns the square change in the length of the spring.

    inline double deltaL(double* pos, const int k) {

        return euclDist(pos, k) - RESTLEN;

    }

    // euclDist is a Euclidean distance calculator that works with the
    // pointer nPtr declared in operator(). It takes into account the periodic
    // boundary condition for the system. This is the function that incorporates
    // strain into the network.
    // 
    // Note that we use jMax + 1 and iMax + 1 in these expressions because they
    // are equivalent to the size of the network netSize. 

    inline double euclDist(const double* pos, const int k) {

        double x = (pos[2 * k] + xshift(k) - pos[0]);
        double y = (pos[2 * k + 1] + yshift(k) - pos[1]);

        return sqrt(x * x + y * y);

    }

    double xshift(const int &k);

    double yshift(const int &k);

};

#endif /*NETWORK_H_*/