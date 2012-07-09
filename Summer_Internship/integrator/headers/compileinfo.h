#ifndef COMPILEINFO_H_
#define COMPILEINFO_H_
// compileinfo.h
// -------------
// compileinfo.h contains information specific to the computer the program is 
// being compiled on. This controls output path information primarily.

#include <string>

const std::string output_path = "/scratch/gpfs/myucht/";

// Rest length for springs.
extern const double RESTLEN = 1;
// Viscosity of the fluid.
extern const double ETA = 0.1;
// Radius for Stokes' drag.
extern const double RADIUS = 0.1;
// Young's modulus for springs.
extern const double YOUNGMOD = 1;

#endif /*COMPILEINFO_H_*/
