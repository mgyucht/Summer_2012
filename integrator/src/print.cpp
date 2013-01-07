// print.cpp
// ---------
// print.cpp contains the class Printer to convert a set of ints, doubles and arrays
// into text files, specifically for integrator.cpp.

#include <boost/lexical_cast.hpp>
#include "print.h"

double Printer::affposx(int r, int c)
{
      return (c + r / 2.0 + affdel / 2.0 * ((2.0 * r) / (netSize - 1.0) - 1.0));
}

double Printer::affposy(int r)
{
      return sqrt(3.0) / 2.0 * r;
}

void Printer::printPos(std::string posFileName) {

    std::ofstream posFile(posFileName.c_str(), std::ios::trunc);

    if (posFile.is_open()) {

        posFile << "NetSize,Strain,YoungMod,pbond,Spr1,Spr2,Spr3" << std::endl;
        posFile << netSize << "," << affdel << "," << YOUNGMOD << ","
            << p << ",1,1,1" << std::endl;

        for (int i = 0; i < netSize; i++) {

            for (int j = 0; j < netSize; j++) {

                // Print row, col, position, sprstiff, rlen information to file.

                posFile << i << "," << j << "," << pos[(i * netSize + j) * 2] - affposx(i, j)
                    << "," << pos[(i * netSize + j) * 2 + 1] - affposy(i);

                for (int k = 0; k < 3; k++)
                    posFile << "," << spr[i][j][k];

                posFile << "\n";
            }

        }

    }

    posFile.close();
}

void Printer::printNonAff(std::string nonaffFileName, int i, double str_rate) {

    if (i == 0)
        clearNonAffFile(nonaffFileName);

    std::ofstream nonaffFile(nonaffFileName.c_str(), std::ios::app);

    if (nonaffFile.is_open()) {
        if (i == 0) {
            nonaffFile << p << "," << netSize << "," << TIMESTEP << std::endl;
        } else {
            double nonaff = nonAffinity(pos);
            double nonaffdd = nonAffinity_dd(pos, del, str_rate);
            nonaffFile << i * TIMESTEP << "," << affdel << ","
                << nonaff << "," << str_rate << "," << nonaffdd << std::endl;
        }
    }

    nonaffFile.close();

}

void Printer::clearNonAffFile(std::string nonaffFileName)
{
  std::ofstream nonaffFile(nonaffFileName.c_str(), std::ios::trunc);
  nonaffFile.close();
}

void Printer::printEnergy(std::string energyFileName, const double &newEnergy) {

    std::ofstream engFile(energyFileName.c_str(), std::ios::app);

    if (engFile.is_open()) {

        engFile << newEnergy << "," << p << "," << affdel << std::endl;

    }

    engFile.close();

}

void Printer::printStress(std::string stressFileName, const double* stress_array,
        const double* strain_array) {

    std::ofstream stressFile(stressFileName.c_str(), std::ios::trunc);

    if (stressFile.is_open()) {

        for (int i = 0; i < n_time_steps; i += fs) {

            std::string time = boost::lexical_cast<std::string>(i * TIMESTEP);

            stressFile << stress_array[i] << "," << strain_array[i] << ",";
            stressFile << time << std::endl;

        }

    }

    stressFile.close();

}
