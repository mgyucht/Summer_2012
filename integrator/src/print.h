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

enum {num_data = 1000};

extern const double YOUNGMOD;
extern int netSize;
extern double strain;
extern double TIMESTEP;
extern int frame_sep;

struct Printer {

    double p;
    double n_time_steps;
    double *pos;
    double *del;
    double ***spr;

    Printer(const Network &net, const double &pp, const double &nts) :
        p(pp),
        n_time_steps(nts),
        pos(net.pos),
        del(net.delta),
        spr(net.spring) {}

    void printPos(std::string /*fileName*/);

    void printNonAff(std::string /*nonaffFileName*/, int /* time */, double /* str_rate */);
    void clearNonAffFile(std::string /*nonaffFileName*/);

    void printEnergy(std::string /*fileName*/, const double & /*newEnergy*/);

    // Prints the stress of the network for every tenth time point.
    void printStress(std::string /*fileName*/, const double* stress, const double* strain);

};

#endif /*PRINT_H_*/
