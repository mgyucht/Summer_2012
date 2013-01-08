/* integrator.cpp
 * --------
 *
 * integrator.cpp is the client file in simulating the spring networks. This
 * program works by integrating the equations of motion over all of the nodes
 * in the network. It also takes into account the viscous force from the
 * shearing of a liquid.
 *
 * Author: Miles Yucht
 * Date: Mon July 07 2012
 */

#include <string>
#include <ctime>
#include <cstdlib>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <boost/program_options.hpp>

namespace po = boost::program_options;

#include "utils.h"
#include "network.h"
#include "nonaffinity.h"
#include "print.h"
#include "options.h"

// Rest length for springs.
const double RESTLEN = 1.0;
// Viscosity of the fluid.
const double ETA = 1.0e0;
// Radius for Stokes' drag.
const double RADIUS = 1.0;
// Young's modulus for springs.
const double YOUNGMOD = 1.0;
// Time step for the simulation.
double TIMESTEP;
double affdel = 0.0;

// NOTE: One unit of time in this simulation is equivalent to 10^-5 seconds in
// reality: 1 s^* = 10^-5 s.

int netSize;          // Network size
double strain;        // Network strain
int frame_sep;     // Number of steps between generated output.

int main (int argc, char *argv[])
{

    Options myOptions;
#ifdef DEBUG
    std::cout << "Calling setup_options" << std::endl;
#endif
    int ret = setup_options(argc, argv, myOptions);
#ifdef DEBUG
    std::cout << "setup_options returned " << ret << std::endl;
#endif
    if (ret)
      return ret;

    int prngseed = myOptions.prngseed,    // Random number generator seed for springs (0)
        nTimeSteps = myOptions.nTimeSteps,  // Number of time steps to simulate (200000)
        steps_per_oscillation,              // Exactly what you think it is
        out_per_oscillation = myOptions.out_per_oscillation, // How many times to output per oscillation
        num_osc = myOptions.num_osc, // Number of oscillations
        motors = myOptions.motors;      // Use motors (1)

    double pBond = myOptions.pBond,             // Bond probability (0.8)
           strRate = myOptions.strRate,           // Strain rate (1.0 Hz*)
           temp = myOptions.temp,              // Temperature of the system
           initStrain = myOptions.initStrain,        // Magnitude of strain (0.01)
           test_step = myOptions.test_step,         // Time step candidate from strain rate
           max_time_step = 0.1; // Maximum time step (0.3 s*) [Constant]

    std::string energyFileName = myOptions.energyFileName, // Energy file name
           posFileName = myOptions.posFileName,    // Position file name
           nonaffFileName = myOptions.nonaffFileName, // Nonaffinity file name
           stressFileName = myOptions.stressFileName, // Stress file name
           output_path = myOptions.output_path,    // Output path for simulation
           config_file = myOptions.config_file,    // Name and location of config file
           job = myOptions.job,            // Job (only used on della) (0)
           extension = ".txt"; // File extension (".txt")

    // Set the time step.

    test_step = 2 * PI / (1000 * strRate);
    TIMESTEP = test_step < max_time_step ? test_step : max_time_step;
    steps_per_oscillation = (int) (strRate > 1e-15 
                                     ? (2 * PI / (strRate * TIMESTEP)) 
                                     : 1000);

    nTimeSteps = steps_per_oscillation * num_osc;
    frame_sep = steps_per_oscillation / out_per_oscillation;

#ifdef DEBUG
    printf("Time step for simulation: %.3g\n"
           "Steps per oscillation:    %.3g\n"
           "Number of time steps:     %d\n", TIMESTEP, steps_per_oscillation, 
           nTimeSteps);
#endif

    // Set the file paths.

#ifdef DELLA3
    std::string root_path = output_path + "/" + job;
#else
    std::string root_path = output_path;
#endif

    bool print_array[4];
    print_array[0] = static_cast<bool>(posFileName.compare(""));
    print_array[1] = static_cast<bool>(nonaffFileName.compare(""));
    print_array[2] = static_cast<bool>(stressFileName.compare(""));
    print_array[3] = static_cast<bool>(energyFileName.compare(""));

    std::string stressFilePath = root_path + "/" + stressFileName + extension;
    std::string energyFilePath = root_path + "/" + energyFileName + extension;
    std::string nonaffFilePath = root_path + "/" + nonaffFileName + extension;

#ifdef DEBUG
    printf("Output path: %s\n"
           "Stress file: %s\n"
           "Energy file: %s\n"
           "Nonaffinity file: %s\n", output_path.c_str(), stressFilePath.c_str(), 
           energyFilePath.c_str(), nonaffFilePath.c_str());
#endif

    // Initialize PRNG.

    if (prngseed == 0)
    {
        unsigned int seed = (unsigned int) time(NULL);
#ifdef DEBUG
        printf("Random seed: %d\n", seed);
#endif
        srand( seed );
    } else
    {
#ifdef DEBUG
        printf("Random seed: %d\n", prngseed);
#endif
        srand( prngseed );
    }

    // Now that those are parsed, we can start to generate our network.

    double *position = new double [2 * netSize * netSize];
    double *delta = new double [2 * netSize * netSize];
    double *stress_array = new double [nTimeSteps];
    double ***sprstiff = new double **[netSize];
    double ***netForces = new double **[netSize];

    for (int i = 0; i < netSize; i++)
    {
        sprstiff[i] = new double *[netSize];
        netForces[i] = new double *[netSize];

        for (int j = 0; j < netSize; j++)
        {
            // x-coordinate
            position[(i * netSize + j) * 2] = RESTLEN * (i / 2.0 + j);

            // y-coordinate
            position[(i * netSize + j) * 2 + 1] = sqrt(3) / 2 * RESTLEN * i;

            sprstiff[i][j] = stiffVecGen(pBond, 3);
            netForces[i][j] = new double[6];
        }
    }

#ifdef DEBUG
    printf("position, delta, stress_array, sprstiff, and netForces are all"
           " allocated.\n");
#endif

    // If the strain magnitude is gamma * network_height, the actual strain on
    // the network is 2 * gamma. Therefore, we halve gamma before making the
    // strain array so that requesting a simulation with a certain strain
    // results in the network with that strain and not double that strain.

    initStrain *= 1 / 2.0;
    double* strain_array = new double[nTimeSteps];
    double* strain_rate = new double[nTimeSteps];

    for (int i = 0; i < nTimeSteps; i++)
    {
        strain_rate[i]  = initStrain * strRate * cos(strRate * i * TIMESTEP);
    }

    Network myNetwork(position, delta, sprstiff, netForces);
    Printer myPrinter(myNetwork, pBond, nTimeSteps, frame_sep);
    Motors myMotors(sprstiff);

#ifdef DEBUG
    printf("myNetwork, myPrinter, myMotors, strain_array, and strain_rate all"
           " allocated.\n");
#endif

    // Integrate motion over the nodes.
    //
    // Technical note: Stress is a rank 2 tensor. However, because our interest
    // is in the shear modulus of the network, we discard the isotropic
    // elements of the stress. Therefore, because this is a 2-D network, there
    // is only one real term of interest, because σ_xy = σ_yx. This is the
    // number that is output by calcStress.

    for (int i = 0; i < nTimeSteps; i++) {
        strain_array[i] = affdel * 2 / (sqrt(3.0) / 2.0 * netSize);;
        // Calculate the net forces in the network.

        if (motors != 0)
            myNetwork.getNetForces(myMotors);
        else
            myNetwork.getNetForces();

        // Calculate the stress of the network.

        stress_array[i] = myNetwork.calcStress(strain_rate[i]);

        // Quit if stress_array[i] is nan.

        if (stress_array[i] != stress_array[i])
        {
            printf("Stress has gone to NaN.\n");
            printf("p = %.2g, w = %.4g, N = %d, e = %.2g\n", pBond, strRate, netSize, initStrain);
            myPrinter.printStress(stressFilePath.c_str(), stress_array, strain_array);
            return 2;
        }

        // If a filename is specified, print the positions of the nodes.

        if (i % frame_sep == 0)
        {
          std::cout << "/" << std::flush;
          if (print_array[0]) { // Position data
            int frame = i / frame_sep;
            std::string iter = boost::lexical_cast<std::string>(frame);
            std::string posFilePath = root_path + posFileName + "_" + iter + extension;
            myPrinter.printPos(posFilePath.c_str());
          }
          if (print_array[1]) { // Time-varying nonaffinity.
            myPrinter.printNonAff(nonaffFilePath.c_str(), i, strain_rate[i > 0 ? i - 1 : 0]);
          }
        }

        // Simulate the movement for this time step.

        myNetwork.moveNodes(strain_rate[i], temp);
    }
    printf("\n");

    // The boolean variables defined above determine whether or not to print
    // this information.

    if (print_array[2])
    {
      myPrinter.printStress(stressFilePath.c_str(), stress_array, strain_array);
    }

    if (print_array[3])
    {
      myPrinter.printEnergy(energyFilePath.c_str(), myNetwork()); // Energy
    }

    // Cleanup
    delete[] position;
    delete[] stress_array;

    for (int i = 0; i < netSize; i++) {
        for (int j = 0; j < netSize; j++) {
            delete[] sprstiff[i][j];
            delete[] netForces[i][j];
        }
        delete[] sprstiff[i];
        delete[] netForces[i];
    }

    delete[] sprstiff;
    delete[] netForces;
    delete[] strain_rate;
    delete[] strain_array;

    return 0;
}
