#ifndef PRINT_H_
#define PRINT_H_

// print.h
// -------
// Header file for print.cpp. Includes function prototypes for printing array
// values and other information to text files.

#include <string>
#include <cstring>
#include <fstream>
#include <sstream>
#include "nonaffinity.h"
#include "network.h"

extern int netSize;
extern double strain;

struct Printer {
    
    double yMod;
    double p;
    double n_time_steps;
    double time_step_size;
    double* pos;
    double*** spr;
    
    Printer(const Network &net, const double &yyMod, const double &pp, 
            const double &nts) : 
        yMod(yyMod), 
        p(pp), 
        n_time_steps(nts) {
        
        this->pos = net.pos;
        this->spr = net.spr;
        time_step_size = net.timestep;
        
    }

    void printPos(std::string /*fileName*/);
    
    void printNonAff(std::string /*fileName*/);
    
    void printEnergy(std::string /*fileName*/, const double & /*newEnergy*/);
    
    void printStress(std::string /*fileName*/, const double* stress, const double* strain);
    
};

#endif /*PRINT_H_*/
