/*
 * printer.cxx
 * -----------
 *
 *  printer.cxx contains the implementation of the Printer class.
 */

#include "printer.hxx"
#include <boost/lexical_cast.hpp>
#include <iostream>

void Printer::print_positions(Simulation sim, int timestep)
{
  pos_file_name = pos_file_name + "_" + boost::lexical_cast<string>(timestep);
  std::vector<Node> *temp_node_vec = &(sim.node_network);
  for (int i = 0; i < temp_node_vec->size(); i++)
  {
    std::vector<double> node_posget_position()
    position_file << i << "," << 
  }
}
