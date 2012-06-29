#ifndef UTILS_H_
#define UTILS_H_

// utils.h
// Simple utilities for energy minimizing code.

#include <stdio.h>
#include <stdlib.h>

extern const double YOUNGMOD;

// Utility function used in this structure or others derived from it.
inline void shft2(double &a, double &b, const double c)
{
    a=b;
    b=c;
}

inline void shft3(double &a, double &b, double &c, const double d)
{
    a=b;
    b=c;
    c=d;
}

inline void mov3(double &a, double &b, double &c, const double d, const double e,
        const double f)
{
    a=d;
    b=e;
    c=f;
}

inline void SWAP(double a, double b) {

    double dum = a;
    a = b;
    b = dum;

}

inline double abs(double a) {

    return a >= 0 ? a : -a;

}

inline double MAX(double a, double b) {

    return a >= b ? a : b;

}

inline double SIGN(double a, double b) {

    return b >= 0 ? (a >= 0 ? a : -a) : (a >= 0 ? -a : a);

}

//usageExit sends an error message about the usage of main and exits the
//program.

inline void usageExit() {
    printf("Usage:");
    // -str -size -p -y -n -step -rate -seed --help
    printf(" integrator [options] \n");
    printf("Options: \n");
    printf("  -str <arg> \t\t Adjust the strain (default: 0.01)\n");
    printf("  -size <arg> \t\t Set the network size (default: 20)\n");
    printf("  -p <arg> \t\t Set the bond probability (default: 0.8)\n");
    printf("  -n <arg> \t\t Set the number of time steps (default: 100)\n");
    printf("  -step <arg> \t\t Set the size of a time step (default: 0.01)\n");
    printf("  -rate <arg> \t\t Set the rate of oscillations (defautl: 1.0)\n");
    printf("  -seed <arg> \t\t Set the seed for random-number generation ");
    printf("(0 for time-based) (default: 0)\n");
    printf("  -energy-fn <arg> \t Set the energy file name (default: ../output/energy_data.txt)\n");
    printf("  -position-fn <arg> \t Set the position file name (default: ../output/position_data.txt)\n");
    printf("  -stress-fn <arg> \t Set the stress file name (default: ../output/stress_data.txt)\n");
    printf("  -nonaff-fn <arg> \t Set the non-affinity file name (default: ../output/nonaff_data.txt)\n");
    printf("  --help \t\t Display this help text\n\n");
    printf("For bug reporting information, contact ME!\n");
    exit(EXIT_FAILURE);
}

//randDouble returns a random number in the range [low, high).
inline double randDouble(int low, int high) {

    double x = low + (high - low) * ((double) rand() / ((double) RAND_MAX + 1));
    return x;

}

//stiffVecGen returns a vector containing three numbers, corresponding to the
//spring constants for the node. These numbers may be either yMod or 0.

inline double *stiffVecGen(double prob, int numSprings) {

    double *ret = new double[numSprings];

    for (int i = 0; i < numSprings; i++)
        ret[i] = randDouble(0, 1) > prob ? 0 : YOUNGMOD;

    return ret;

}

#endif /*UTILS_H_*/
