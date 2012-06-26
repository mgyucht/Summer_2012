// print.cpp
// ---------
// print.cpp contains the class Printer to convert a set of ints, doubles and arrays
// into text files, specifically for integrator.cpp.

#include <string>
#include <cstring>
#include <fstream>
#include "nonaffinity.h"
#include "print.h"

void Printer::printPos(std::string posFileName) {

    std::ofstream posFile(posFileName.c_str(), std::ios::trunc);

    if (posFile.is_open()) {
 
        posFile << "Network Size,Strain,Young's Modulus,Probability" << std::endl;
        posFile << netSize << "," << strain << "," << yMod << ","
            << p << std::endl;

        for (int i = 0; i < netSize; i++) {

            for (int j = 0; j < netSize; j++) {

                // Print row, col, position, sprstiff, rlen information to file.

                posFile << i << "," << j << "," << pos[(i * netSize + j) * 2]
                    << "," << pos[(i * netSize + j) * 2 + 1];

                for (int k = 0; k < 3; k++) 
                    posFile << "," << spr[i][j][k];

                posFile << std::endl;
            }

        }

    }

    posFile.close();
}

void Printer::printNonAff(std::string nonaffFileName) {

    std::ofstream nonaffFile(nonaffFileName.c_str(), std::ios::app);

    if (nonaffFile.is_open()) {
        
        double nonaff = nonAffinity(nonaffFileName.c_str());
        nonaffFile << netSize << "," << strain << "," << p << "," 
            << nonaff << std::endl;
    
    }
    
    nonaffFile.close();
    
}

void Printer::printEnergy(std::string energyFileName, const double &newEnergy) {

    std::ofstream engFile(energyFileName.c_str(), std::ios::app);

    if (engFile.is_open()) {

        engFile << newEnergy << "," << p << "," << strain << std::endl;

    }

    engFile.close();
    
}

void Printer::printStress(std::string stressFileName, const double *stress, 
        const int &steps) {
    
    std::ofstream stressFile(stressFileName.c_str(), std::ios::trunc);

    if (stressFile.is_open()) {
        
        for (int i = 0; i < steps; i++) 

            stressFile << stress[i] << "," << i << std::endl;
    
    }

}
