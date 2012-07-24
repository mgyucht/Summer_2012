// print.cpp
// ---------
// print.cpp contains the class Printer to convert a set of ints, doubles and arrays
// into text files, specifically for integrator.cpp.

#include "print.h"

void Printer::printPos(std::string posFileName) {

    std::ofstream posFile(posFileName.c_str(), std::ios::trunc);

    if (posFile.is_open()) {
 
        posFile << "NetSize,Strain,YoungMod,pbond,Spr1,Spr2,Spr3" << std::endl;
        posFile << netSize << "," << strain << "," << YOUNGMOD << ","
            << p << ",1,1,1" << std::endl;

        for (int i = 0; i < netSize; i++) {

            for (int j = 0; j < netSize; j++) {

                // Print row, col, position, sprstiff, rlen information to file.

                posFile << i << "," << j << "," << pos[(i * netSize + j) * 2]
                    << "," << pos[(i * netSize + j) * 2 + 1];

                for (int k = 0; k < 3; k++) 
                    posFile << "," << spr[i][j][k];

                posFile << "\n";
            }

        }

    }

    posFile.close();
}

void Printer::printNonAff(std::string nonaffFileName, int i) {

    std::ofstream nonaffFile(nonaffFileName.c_str(), std::ios::trunc);

    if (nonaffFile.is_open()) {
        
        double nonaff = nonAffinity(pos, netSize, strain);
        nonaffFile << p << "," << i * TIMESTEP << "," << strain << "," 
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

void Printer::printStress(std::string stressFileName, const double* stress_array, 
        const double* strain_array) {
    
    std::ofstream stressFile(stressFileName.c_str(), std::ios::trunc);

    if (stressFile.is_open()) {
        
        for (int i = 0; i < n_time_steps; i += 10) {

            std::ostringstream int_to_str(std::ostringstream::out);
            int_to_str << i * TIMESTEP;
            
            stressFile << stress_array[i] << "," << strain_array[i] << ",";
            stressFile << int_to_str.str() << std::endl;
    
        }
        
    }

}
