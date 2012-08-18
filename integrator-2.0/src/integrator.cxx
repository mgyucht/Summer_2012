/*
 * integrator.cxx
 * --------------
 *
 *  integrator.cxx contains the source of the main program to run the spring
 *  simulation.
 */

#include "printer.hxx"
#include <cmath>
#include <boost/program_options.hpp>

#define PI (3.14159265358979)

namespace po = boost::program_options;

void Simulation::initialize_strain()
{
  double sf = Simulation::strain_frequency;
  double sa = Simulation::strain_amplitude;
  double test_step = 2 * PI / (1000 * Simulation::strain_frequency);
  Simulation::dt = test_step < max_time_step ? test_step : max_time_step;
  Simulation::steps_per_oscillation =
    (int) (Simulation::strain_frequency > 1e-15 ? 2 * PI
             / (Simulation::strain_frequency * Simulation::dt)
             : 200000);
  double total_steps = Simulation::steps_per_oscillation
                        * Simulation::num_of_oscillations;
  for (int i = 0; i < total_steps; i++)
  {
    strain[i] = sa * sin(sf * i * Simulation::dt);
    strain_rate[i] = sa * sf * cos(sf * i * Simulation::dt);
  }
}

int get_simulation_info(Simulation sim, Printer prt, int argc, char* argv[])
{
  std::string config_file;
  std::string job;
  int netSize;

  po::options_description general("General options");
  general.add_options()
    ("help,h", "show this help text")
    ("config,c", po::value<std::string>(&config_file)->default_value(
                                     ".integratorconf"), "set the config file")
    ("netsize,z", po::value<int>(&netSize)->default_value(20),
     "set network dimensions")
    ("probability,p", po::value<double>(&(sim.p_bond))->default_value(0.8),
     "set bond probability")
    ("rate,r", po::value<double>(&(sim.strain_frequency))->default_value(1.0),
     "set oscillation frequency")
    ("strain,e", po::value<double>(&(sim.strain_amplitude))->default_value(0.01),
     "set initial strain")
    ("temp,t", po::value<double>(&(sim.temp))->default_value(0.0),
     "set temperature")
    ("prng", po::value<int>(&(sim.prng_seed))->default_value(0),
     "set PRNG seed")
    ("job,j", po::value<std::string>(&job)->default_value(""), "set job directory")
    ("motors,m", po::value<bool>(&(sim.use_motors))->default_value(false), "enable motors")
    ;

  po::options_description filename("Filename options");
  filename.add_options()
    ("aff-fn", po::value<std::string>(&(prt.nonaff_file_name))->default_value(""),
     "set non-affinity data file name")
    ("pos-fn", po::value<std::string>(&(prt.pos_file_name))->default_value(""),
     "set position data file name")
    ("st-fn", po::value<std::string>(&(prt.stress_file_name))->default_value(""),
     "set stress data file name")
    ;

  po::options_description config("Configuration");
  config.add_options()
    ("output", po::value<std::string>(&(prt.output_path))->default_value(""),
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
    std::cout << cmdline_options << std::endl;
    return 10;
  }

  // Read parameters from the config file.

  std::ifstream ifs(config_file.c_str());

  if (!ifs)
  {
    std::cout << "Couldn't open config file. Please create one and store the"
      << " output information in it." << std::endl;
    return 1;
  } else
  {
    store(parse_config_file(ifs, config_file_options), vm);
    notify(vm);
  }

  // Requires an output_path to be assigned in config file.

  if (prt.output_path.empty())
  {
    std::cout << "You must specify an output path in the config file."
      << std::endl;
    return 1;
  }

}

int main(int argc, char* argv[])
{
  Simulation simulation;
  Printer printer;
  int ret = get_simulation_info(simulation, printer, argc, argv);
  if (ret != 0)
  {
    return ret;
  }
  simulation.initialize_strain();

  return 0;
}
