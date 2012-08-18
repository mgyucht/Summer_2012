/*
 * integrator.hxx
 * --------------
 *
 *  integrator.hxx defines the Simulation class, which incorporates all of the
 *  parameters of the simulation.
 */

#include "node.hxx"
#include <vector>

const double max_time_step = 0.1;

class Simulation
{
  public:
    /* Default initializer */
    Simulation()
    {
      prng_seed = 0;
      p_bond = 0.8;
      eta = 1.0;
      youngs_mod = 1.0;
      temp = 0.0;
      strain_amplitude = 0.1;
      strain_frequency = 1.0;
      initialize_strain();
    }

    void initialize_strain();

    bool use_motors;
    int prng_seed;
    int steps_per_oscillation;
    int out_per_oscillation;
    int num_of_oscillations;
    double eta;
    double youngs_mod;
    double p_bond;
    double strain_frequency;
    double strain_amplitude;
    double temp;
    double dt;
    std::vector<double> strain;
    std::vector<double> strain_rate;
    std::vector<double> stress;
    std::vector<double> network_dimensions;
    std::vector<Node> node_network;
};
