#ifndef NETWORK_H_
#define NETWORK_H_

// network.h
// -------
//
// network.h is the header for network.cpp, which contains the essential methods 
// for integrate.cpp. It defines the struct Network.

#define xadj (pos[2 * k] + xshift)
#define yadj (pos[2 * k + 1] + yshift)
 
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

    Network(double* ppos, double*** espr, double*** vvels, double*** fforces,
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

};

#endif /*NETWORK_H_*/
