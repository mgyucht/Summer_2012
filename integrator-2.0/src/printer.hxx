/*
 * printer.hxx
 * -----------
 *
 *  printer.hxx provides the Printer class.
 */

#include "integrator.hxx"
#include <string>
#include <vector>

class Printer
{
  public:
    Printer()
    {
      job = "0";
      extension = ".txt";
      stress_file_name = "stress_data";
      nonaff_file_name = "nonaff_data";
      pos_file_name = "pos_data";
    }

    void print_positions(Simulation /* sim */, int /* timestep*/);
    void clear_data_file();
    void print_data_file(Simulation /* sim */);

    std::string job;
    std::string extension;
    std::string stress_file_name;
    std::string nonaff_file_name;
    std::string pos_file_name;
    std::string output_path;
};
