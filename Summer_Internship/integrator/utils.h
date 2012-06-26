#ifndef UTILS_H_
#define UTILS_H_

// utils.h
// Simple utilities for energy minimizing code.

#include <stdio.h>

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

void usageExit() {
    printf("usage:");
    printf("\n   integrator [-size <network size>] [-str <initial strain>] [-rate <strain rate>]\n");
    printf("              [-p <bond probability>] [-y <young's modulus for spring>] \n");
    printf("              [-seed <PRNG seed>] [-step <time step>]\n");
    exit(EXIT_FAILURE);
}

//randDouble returns a random number in the range [low, high).
double randDouble(int low, int high) {

    double x = low + (high - low) * ((double) rand() / ((double) RAND_MAX + 1));
    return x;

}

//stiffVecGen returns a vector containing three numbers, corresponding to the
//spring constants for the node. These numbers may be either yMod or 0.

double *stiffVecGen(double prob, double yMod, int numSprings) {

    double *ret = new double[numSprings];

    for (int i = 0; i < numSprings; i++)
        ret[i] = randDouble(0, 1) > prob ? 0 : yMod;

    return ret;

}

#endif /*UTILS_H_*/
