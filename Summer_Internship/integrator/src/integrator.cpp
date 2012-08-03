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

using namespace std;

// Rest length for springs.
const double RESTLEN = 1.0;
// Viscosity of the fluid.
const double ETA = 1.0e0;
// Radius for Stokes' drag.
const double RADIUS = 0.1;
// Young's modulus for springs.
const double YOUNGMOD = 1.0;
// Time step for the simulation.
double TIMESTEP;

// NOTE: One unit of time in this simulation is equivalent to 10^-5 seconds in
// reality: 1 s^* = 10^-5 s.

int netSize;          // Network size
double strain;        // Network strain
int frame_sep;     // Number of steps between generated output.

int main (int argc, char *argv[])
{
    int prngseed,    // Random number generator seed for springs (0)
        nTimeSteps,  // Number of time steps to simulate (200000)
        steps_per_oscillation, // Exactly what you think it is
        out_per_oscillation = 32, // How many times to output per oscillation
        motors;      // Use motors (1)

    double pBond,             // Bond probability (0.8)
           strRate,           // Strain rate (1.0 Hz*)
           temp,              // Temperature of the system
           initStrain,        // Magnitude of strain (0.01)
           test_step,         // Time step candidate from strain rate
           max_time_step = 0.3; // Maximum time step (0.3 s*) [Constant]

    string energyFileName, // Energy file name
           posFileName,    // Position file name
           nonaffFileName, // Nonaffinity file name
           stressFileName, // Stress file name
           output_path,    // Output path for simulation
           config_file,    // Name and location of config file
           job,            // Job (only used on della) (0)
           extension = ".txt"; // File extension

    // Set up command-line parameters and config information.

    po::options_description general("General options");
    general.add_options()
        ("help,h", "show this help text")
        ("config,c", po::value<string>(&config_file)->default_value(
                            ".integratorconf"), "set the config file")
        ("netsize,z", po::value<int>(&netSize)->default_value(20),
             "set network dimensions")
        ("probability,p", po::value<double>(&pBond)->default_value(0.8),
             "set bond probability")
        ("rate,r", po::value<double>(&strRate)->default_value(1.0),
             "set oscillation frequency")
        ("strain,e", po::value<double>(&initStrain)->default_value(0.01),
             "set initial strain")
        ("temp,t", po::value<double>(&temp)->default_value(0.0),
             "set temperature")
        ("prng", po::value<int>(&prngseed)->default_value(0),
             "set PRNG seed")
        ("job,j", po::value<string>(&job)->default_value(""), "set job directory")
        ("motors,m", po::value<int>(&motors)->default_value(0), "enable motors")
        ;

    po::options_description filename("Filename options");
    filename.add_options()
        ("en-fn", po::value<string>(&energyFileName)->default_value(""),
             "set energy data file name")
        ("aff-fn", po::value<string>(&nonaffFileName)->default_value(""),
             "set non-affinity data file name")
        ("pos-fn", po::value<string>(&posFileName)->default_value(""),
             "set position data file name")
        ("st-fn", po::value<string>(&stressFileName)->default_value(""),
             "set stress data file name")
        ;

    po::options_description config("Configuration");
    config.add_options()
        ("output", po::value<string>(&output_path)->default_value(""),
                             "set output path")
        ;

    po::options_description cmdline_options;
    cmdline_options.add(general).add(filename);

    po::options_description config_file_options;
    config_file_options.add(filename).add(config);

    // Make the variables_map object.
    po::variables_map vm;

    // Read parameters from the command line.

    store(po::parse_command_line(argc, argv, cmdline_options), vm);
    notify(vm);

    if (vm.count("help"))
    {
        cout << cmdline_options << endl;
        return 0;
    }

    // Read parameters from the config file.

    ifstream ifs(config_file.c_str());

    if (!ifs)
    {
        cout << "Couldn't open config file. Please create one and store the"
            << " output information in it.\n";
        return 1;
    } else
    {
        store(parse_config_file(ifs, config_file_options), vm);
        notify(vm);
    }

    // Requires an output_path to be assigned in config file.

    if (output_path.empty())
    {
        cout << "You must specify an output path in the config file.\n";
        return 1;
    }

    // Set the time step.

    test_step = 2 * PI / (1000 * strRate);
    TIMESTEP = test_step < max_time_step ? test_step : max_time_step;
    steps_per_oscillation = (int) (strRate > 1e-15 ? (2 * PI / (strRate * TIMESTEP)) : 1000);

    nTimeSteps = steps_per_oscillation * 5;
    frame_sep = steps_per_oscillation / out_per_oscillation;

    // Set the file paths.

    string root_path = output_path + job + "/";

    bool print_array[4];
    print_array[0] = static_cast<bool>(posFileName.compare(""));
    print_array[1] = static_cast<bool>(nonaffFileName.compare(""));
    print_array[2] = static_cast<bool>(stressFileName.compare(""));
    print_array[3] = static_cast<bool>(energyFileName.compare(""));

    string stressFilePath = root_path + stressFileName + extension;
    string energyFilePath = root_path + energyFileName + extension;

    // Initialize PRNG.

    if (prngseed == 0)
        srand( time(NULL) );
    else
        srand( prngseed );

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

    // If the strain magnitude is gamma * network_height, the actual strain on
    // the network is 2 * gamma. Therefore, we halve gamma before making the
    // strain array so that requesting a simulation with a certain strain
    // results in the network with that strain and not double that strain.

    initStrain /= 2;
    double* strain_array = new double[nTimeSteps];
    double* strain_rate = new double[nTimeSteps];

    for (int i = 0; i < nTimeSteps; i++)
    {
        strain_array[i] = initStrain * sin(strRate * i * TIMESTEP);
        strain_rate[i]  = initStrain * strRate * cos(strRate * i * TIMESTEP);
    }

    Network myNetwork(position, delta, sprstiff, netForces);
    Printer myPrinter(myNetwork, pBond, nTimeSteps);
    Motors myMotors(sprstiff);

    {
      string nonaffFilePath = root_path + nonaffFileName + extension;
      myPrinter.clearNonAffFile(nonaffFilePath.c_str());
    }

    // Integrate motion over the nodes.
    //
    // Technical note: Stress is a rank 2 tensor. However, because our interest
    // is in the shear modulus of the network, we discard the isotropic
    // elements of the stress. Therefore, because this is a 2-D network, there
    // is only one real term of interest, because σ_xy = σ_yx. This is the
    // number that is output by calcStress.

    for (int i = 0; i < nTimeSteps; i++)
    {
        strain = strain_array[i];

        // Calculate the net forces in the network.

        if (motors != 0)
            myNetwork.getNetForces(myMotors);
        else
            myNetwork.getNetForces();

        // Calculate the stress of the network.

        stress_array[i] = myNetwork.calcStress(strain_rate[i]);

        // Test
        /*
        if (i % frame_sep == 0)
        {
            printf("%g\n", stress_array[i]);
        }
        */

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
          if (print_array[0]) // Position data
          {
            string iter = boost::lexical_cast<string>(i);
            string posFilePath   = root_path + posFileName + "_" + iter + extension;
            myPrinter.printPos(posFilePath.c_str());
          }
          // cout << "/" << flush;
          printf("%g, %g\n", nonAffinity(position), nonAffinity_dd(position, delta, strain_rate[i]));
          if (print_array[1]) // Time-varying nonaffinity.
          {
            string nonaffFilePath = root_path + nonaffFileName + extension;
            myPrinter.printNonAff(nonaffFilePath.c_str(), i, strain_rate[i]);
          }
        }

        // Simulate the movement for this time step.

        myNetwork.moveNodes(strain_rate[i], temp);
    }
    printf("\n");

    // The boolean variables defined above determine whether or not to print
    // this information.

    if (print_array[2]) myPrinter.printStress(stressFilePath.c_str(), stress_array,
            strain_array);

    if (print_array[3]) myPrinter.printEnergy(energyFilePath.c_str(), myNetwork()); // Energy

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
