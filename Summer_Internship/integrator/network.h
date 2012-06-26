#ifndef NETWORK_H_
#define NETWORK_H_

// network.h
// -------
//
// network.h is the header for network.cpp, which contains the essential methods 
// for integrate.cpp. It defines the struct Network.

#include <math.h>
 
extern const double RESTLEN;

extern int netSize;
extern double strain;

const double DEL = 1E-11;
const double MASS = 1;

struct Network {
    
    double* pos;
    double*** spr;
    double*** vels;
    double*** forces;
    double timestep;
    
    int iMax, jMax;
    
    bool isiMax, isjMax, isiMin, isjMin;

    Network(double* ppos, double*** sspr, double*** vvels, double*** fforces,
        double ttimestep) : 
        pos(ppos), 
        spr(sspr), 
        vels(vvels), 
        forces(fforces),
        timestep(ttimestep) {
            
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
    
    void getNetForces();
    
    double calcStress();
    
    void moveNodes();
    
    private:
    
    // deltaLSqrd returns the square change in the length of the spring.

    inline double deltaLSqrd(double* pos, const int row, const int col, const int k) {

        double lijk = euclDist(pos, k);

        return (lijk - RESTLEN) * (lijk - RESTLEN);

    }

    double euclDist(const double* /* pos */, const int /* k */);
    
    double arctan2(double* pos, const int k);

};

#endif /*NETWORK_H_*/
