#ifndef OPTIONS_H
#define OPTIONS_H
/*
 * options.h
 * ---------
 *  options.h defines a structure allowing me to encapsulate the retrieval of 
 *  options in the .integratorconf file, so integrator.cpp doesn't get too 
 *  cluttered.
 */

#include <string>

class Options
{
  public:

    int prngseed,    // Random number generator seed for springs (0)
        nTimeSteps,  // Number of time steps to simulate (200000)
        out_per_oscillation, // How many times to output per oscillation
        num_osc, // Number of oscillations
        motors;      // Use motors (1)

    double pBond,             // Bond probability (0.8)
           strRate,           // Strain rate (1.0 Hz*)
           temp,              // Temperature of the system
           initStrain,        // Magnitude of strain (0.01)
           test_step; // Maximum time step (0.3 s*) [Constant]

    std::string energyFileName, // Energy file name
           posFileName,    // Position file name
           nonaffFileName, // Nonaffinity file name
           stressFileName, // Stress file name
           output_path,    // Output path for simulation
           config_file,    // Name and location of config file
           job,            // Job (only used on della) (0)
           extension; // File extension
  
};

int setup_options(int argc, char *argv[], Options &myOptions);

#endif /* OPTIONS_H */
