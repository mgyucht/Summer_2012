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
#include <time.h>
#include <stdlib.h>
#include <math.h>
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
const double ETA = 0.1;
// Radius for Stokes' drag.
const double RADIUS = 0.1;
// Young's modulus for springs.
const double YOUNGMOD = 1.0;
// Time step for the simulation.
double TIMESTEP;

// NOTE: One unit of time in this simulation is equivalent to 10^-5 seconds in
// reality: 1 s^* = 10^-5 s.

int netSize;
double strain;

int main (int argc, char *argv[]) {

    // Initialization and default values. The `job' variable is specifically
    // for della.
    int prngseed, nTimeSteps, job;
    double pBond, strRate, currentTime = 0.0, initStrain;
    string energyFileName, posFileName, nonaffFileName, stressFileName,
           output_path, config_file, extension = ".txt";
     
    po::options_description general("General options");
    general.add_options() 
        ("help,h", "show this help text")
        ("config,c", po::value<string>(&config_file)->default_value(".integratorconf"),
             "set the config file")
        ("netsize,z", po::value<int>(&netSize)->default_value(20),
             "set network dimensions")
        ("time,n", po::value<int>(&nTimeSteps)->default_value(20000), 
             "set number of time steps")
        ("probability,p", po::value<double>(&pBond)->default_value(0.8), 
             "set bond probability")
        ("rate,r", po::value<double>(&strRate)->default_value(1.0), 
             "set oscillation frequency")
        ("strain,e", po::value<double>(&initStrain)->default_value(0.01), 
             "set initial strain")
        ("prng", po::value<int>(&prngseed)->default_value(0), 
             "set PRNG seed")
        ("job,j", po::value<int>(&job)->default_value(0), "set job number")
        ;
    
    po::options_description filename("Filename options");
    filename.add_options()
        ("en-fn", po::value<string>(&energyFileName)->default_value("energy_data"),
             "set energy data file name")
        ("aff-fn", po::value<string>(&nonaffFileName)->default_value("nonaff_data"),
             "set non-affinity data file name")
        ("pos-fn", po::value<string>(&posFileName)->default_value("position_data"),
             "set position data file name")
        ("st-fn", po::value<string>(&stressFileName)->default_value("stress_data"),
             "set stress data file name")
        ;
    
    po::options_description config("Configuration");
    config.add_options()
        ("output", po::value<string>(&output_path)->default_value(""), "set output path")
        ;
    
    po::options_description cmdline_options;
    cmdline_options.add(general).add(filename);
    
    po::options_description config_file_options;
    config_file_options.add(filename).add(config);
    
    po::variables_map vm;
    store(po::parse_command_line(argc, argv, cmdline_options), vm);
    notify(vm);
    
    if (vm.count("help")) {
        
        cout << cmdline_options << endl;
        return 0;
    
    }
    
    ifstream ifs(config_file.c_str());
    
    if (!ifs) {
        
        cout << "Couldn't open config file. Please create one and store the" 
            << " output information in it.\n";
        return 1;
        
    } else {
        
        store(parse_config_file(ifs, config_file_options), vm);
        notify(vm);
        
    }
    
    if (output_path.empty()) {
        
        cout << "You must specify an output path in .integratorconf.\n";
        return 1;
        
    }
    
    // Set the timestep. This needs to change for slow strain frequencies
    // (omega < 0.01)
    TIMESTEP = 1 / (1000 * strRate);
    
    string root_path;
    
    if (job > 0) {
        
        // Get filenames ready.
        ostringstream convert;
        convert << job << "/";
        
        root_path = output_path + convert.str();
        
    } else {
        
        root_path = output_path; 
        
    }
    
    // Make the directory if it doesn't exist.
    
    string posFilePath    = root_path + posFileName + extension;
    string nonaffFilePath = root_path + nonaffFileName + extension;
    string stressFilePath = root_path + stressFileName + extension;
    string energyFilePath = root_path + energyFileName + extension;
    
    char* posFileFull    = (char *) (posFilePath).c_str();
    char* nonaffFileFull = (char *) (nonaffFilePath).c_str();
    char* stressFileFull = (char *) (stressFilePath).c_str();
    char* energyFileFull = (char *) (energyFilePath).c_str();
    
    // Initialize PRNG.
    if (prngseed == 0) 

        srand( time(NULL) );
    
    else 

        srand( prngseed );

    // Now that those are parsed, we can start to generate our network.
    double* position = new double [2 * netSize * netSize];
    double* stress_array = new double [nTimeSteps];
    double*** sprstiff = new double** [netSize];
    double*** velocities = new double** [netSize];
    double*** netForces = new double** [netSize];

    for (int i = 0; i < netSize; i++) {

        sprstiff[i] = new double* [netSize];
        velocities[i] = new double* [netSize];
        netForces[i] = new double* [netSize];

        for (int j = 0; j < netSize; j++) {

            // x-coordinate
            position[(i * netSize + j) * 2] = RESTLEN * (i / 2.0 + j);

            // y-coordinate
            position[(i * netSize + j) * 2 + 1] = sqrt(3) / 2 * RESTLEN * i;

            sprstiff[i][j] = stiffVecGen(pBond, 3);
            velocities[i][j] = new double[2];
            netForces[i][j] = new double[6];
            
        }
    }
    
    double* strain_array = new double[nTimeSteps];
    double* strain_rate = new double[nTimeSteps];
    
    for (int i = 0; i < nTimeSteps; i++) {
    
        strain_array[i] = initStrain * sin(strRate * i * TIMESTEP);
        strain_rate[i]   = initStrain * strRate * TIMESTEP * cos(strRate * i * TIMESTEP);
    
    }

    Network myNetwork(position, sprstiff, velocities, netForces);
    Printer myPrinter(myNetwork, pBond, nTimeSteps);
    
    // Integrate motion over the nodes.
    // 
    // Technical note: Stress is a rank 2 tensor. However, because our interest
    // is in the shear modulus of the network, we discard the isotropic
    // elements of the stress. Therefore, because this is a 2-D network, there 
    // is only one real term of interest, because σ_xy = σ_yx. This is the 
    // number that is output by calcStress.
    
    for (int i = 0; i < nTimeSteps; i++) {
    
        currentTime += TIMESTEP;
        strain = strain_array[i];
        
        // Calculate the net forces in the network.
        myNetwork.getNetForces();
        // Calculate the stress of the network.
        stress_array[i] = myNetwork.calcStress();
        // Simulate the movement for this time step.
        myNetwork.moveNodes(strain_rate[i]);
        
    }
    
    // Calculate the energy of the network. Not important for dynamic simulation.
    // double newEnergy = myNetwork();
    
    // Uncomment whichever data you want to have printed out.
    
    // myPrinter.printPos(posFileFull);             // Position
    // myPrinter.printNonAff(nonaffFileFull);       // Non-affinity
    myPrinter.printStress(stressFileFull, stress_array, strain_array);
    // myPrinter.printEnergy(energyFileFull, newEnergy);

    // Cleanup
    delete[] position;
    delete[] stress_array;

    for (int i = 0; i < netSize; i++) {
        for (int j = 0; j < netSize; j++) {
            delete[] sprstiff[i][j];
            delete[] velocities[i][j];
            delete[] netForces[i][j];
        }
        delete[] sprstiff[i];
        delete[] velocities[i];
        delete[] netForces[i];
    }

    delete[] sprstiff;
    delete[] velocities;
    delete[] netForces;
    delete[] strain_rate;
    delete[] strain_array;

    return 0;
}
