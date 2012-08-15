#ifndef _DEBUG_H_
#define _DEBUG_H_

// debug.h
// -------
//
// Defines debugging macros for the preprocessor.

// Set to 1 to enable info for the funcd() energy function.
#define DEBUGENERGY 0
#define DEBUGDIST 0

// Set to 1 to enable info for the funcd.df() gradient function.
#define DEBUGGRADIENT 0

// Set to 1 to enable info for the Frprmn::minimize(double *pp) function.
#define DEBUGFRPRMN 0

// The following macros enable different output from the program.
// PRINTENERGY prints the energy of the system to energy_data.txt.
#define PRINTENERGY 0

// PRINTNONAFFINITY can be enabled only if PRINTPOS is as well.
#define PRINTPOS 1
#define PRINTNONAFFINITY 1

#endif /*_DEBUG_H_*/
