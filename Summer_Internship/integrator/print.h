#ifndef PRINT_H_
#define PRINT_H_
// print.h
// -------
// Header file for print.cpp. Includes function prototypes for printing array
// values and other information to text files.

#include "network.h"

extern int netSize;
extern double strain;

struct Printer {
    
    double yMod;
    double p;
    double* pos;
    double*** spr;
    
    Printer(const Network &net, const double &yyMod, const double &pp) : 
        yMod(yyMod), 
        p(pp) {
        
        this->pos = net.pos;
        this->spr = net.spr;
        
    }

    void printPos(std::string /*fileName*/);
    
    void printNonAff(std::string /*fileName*/);
    
    void printEnergy(std::string /*fileName*/, const double & /*newEnergy*/);
    
    void printStress(std::string /*fileName*/, const double* stress, const int &steps);
    
};

#endif /*PRINT_H_*/
