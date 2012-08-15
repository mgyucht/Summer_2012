#ifndef UTILS_H_
#define UTILS_H_

// utils.h
// Simple utilities for integrator code.

#include <stdlib.h>

extern const double YOUNGMOD;

//randDouble returns a random number in the range (low, high).
inline double randDouble(int low, int high) {

    double x = low + (high - low) * ((double) (rand() + 1) / ((double) RAND_MAX + 1));
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
