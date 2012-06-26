#ifndef DEBUG_H_
#define DEBUG_H_

// debug.h
// -------
//
// Defines debugging macros for the preprocessor.

// The following macros enable different output from the program.
// PRINTENERGY prints the energy of the system to energy_data.txt.
#define PRINTENERGY 0

// PRINTNONAFFINITY can be enabled only if PRINTPOS is as well.
#define PRINTPOS 1
#define PRINTNONAFFINITY 1

#endif /*DEBUG_H_*/
