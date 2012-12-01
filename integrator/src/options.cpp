/*
 * options.cpp
 * ---------
 * Implements options.h functions. Retrieves settings from 
 * .integratorconf.
 */

#include <boost/program_options.hpp>
#include <iostream>
#include <fstream>
#include <stdlib.h>

#include "options.h"

extern int netSize;

int setup_options(int argc, char *argv[], Options &myOpts)
{
    // Set up command-line parameters and config information.
    double *pBond = &(myOpts.pBond);
    double *strRate = &(myOpts.strRate);
    double *initStrain = &(myOpts.initStrain);
    double *temp = &(myOpts.temp);
    int *prngseed = &(myOpts.prngseed);
    int *numosc = &(myOpts.num_osc);
    int *out_per_oscillation = &(myOpts.out_per_oscillation);
    std::string *job = &(myOpts.job);
    int *motors = &(myOpts.motors);
    std::string *energyFileName = &(myOpts.energyFileName);
    std::string *nonaffFileName = &(myOpts.nonaffFileName);
    std::string *stressFileName = &(myOpts.stressFileName);
    std::string *posFileName = &(myOpts.posFileName);
    std::string config_file;
    std::string output_path;
    std::string home = getenv("HOME");

    boost::program_options::options_description general("General options");
    general.add_options()
        ("help,h", "show this help text")
        ("config,c", boost::program_options::value<std::string>(&config_file)->default_value(
                           home + "/.integratorconf"), "set the config file")
        ("netsize,z", boost::program_options::value<int>(&netSize)->default_value(20),
             "set network dimensions")
        ("probability,p", boost::program_options::value<double>(pBond)->default_value(0.8),
             "set bond probability")
        ("rate,r", boost::program_options::value<double>(strRate)->default_value(1.0),
             "set oscillation frequency")
        ("strain,e", boost::program_options::value<double>(initStrain)->default_value(0.01),
             "set initial strain")
        ("temp,t", boost::program_options::value<double>(temp)->default_value(0.0),
             "set temperature")
        ("prng", boost::program_options::value<int>(prngseed)->default_value(0),
             "set PRNG seed")
        ("job,j", boost::program_options::value<std::string>(job)->default_value("0"), "set job directory")
        ("motors,m", boost::program_options::value<int>(motors)->default_value(0), "enable motors")
        ;

    boost::program_options::options_description filename("Filename options");
    filename.add_options()
        ("en-fn", boost::program_options::value<std::string>(energyFileName)->default_value(""),
             "set energy data file name")
        ("aff-fn", boost::program_options::value<std::string>(nonaffFileName)->default_value(""),
             "set non-affinity data file name")
        ("pos-fn", boost::program_options::value<std::string>(posFileName)->default_value(""),
             "set position data file name")
        ("st-fn", boost::program_options::value<std::string>(stressFileName)->default_value(""),
             "set stress data file name")
        ("num-osc", boost::program_options::value<int>(numosc)->default_value(6),
             "set the number of full oscillations")
        ("out-per-osc", boost::program_options::value<int>(out_per_oscillation)->default_value(20),
             "set the number of data points to output per oscillation")
        ;

    boost::program_options::options_description config("Configuration");
    config.add_options()
        ("output", boost::program_options::value<std::string>(&output_path)->default_value(""),
                             "set output path")
        ;

    boost::program_options::options_description cmdline_options;
    cmdline_options.add(general).add(filename);

    boost::program_options::options_description config_file_options;
    config_file_options.add(filename).add(config);

    // Make the variables_map object.
    boost::program_options::variables_map vm;

    // Read parameters from the command line.

    store(boost::program_options::parse_command_line(argc, argv, cmdline_options), vm);
    notify(vm);

    if (vm.count("help"))
    {
      std::cout << cmdline_options << std::endl;
      return 1;
    }

    // Read parameters from the config file.

    std::ifstream ifs(config_file.c_str());

    if (!ifs)
    {
      std::cout << "Couldn't open config file. Please create one and store the"
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
      std::cout << "You must specify an output path in the config file.\n";
      return 1;
    }

    myOpts.config_file = config_file;
    myOpts.output_path = output_path;

    return 0;
}
