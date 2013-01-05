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

    double affdel;
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

        affdel = 0.0;
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


    double xshift(const int &k);

    double yshift(const int &k);

};

#endif /*NETWORK_H_*/
