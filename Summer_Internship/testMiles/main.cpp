/* main.cpp
 * --------
 *
 * main.cpp is the client file in simulating the spring networks. 
 * Usage: main -str <strain> -size <network size> -p <bond probability> -y \
 *            <young's modulus for springs>
 *  
 * Any other usage will result in the program exiting.
 * Author: Miles Yucht
 * Date: Mon June 11 2012
 */

/*
 * Where are things defined? :
 *  RESTLEN -> funcd.h
 */

/*
 *  Working revisions:
 *
 *  Currently, doesn't work.
 *  Only two dimensions.
 *
 */

#include <iostream>
#include <vector>
#include <string>
#include <ctime>
#include <time.h>
#include "nr3.h"
#include "funcd.h"
#include <math.h>
#include "node.h"

using namespace std;

void usageExit();
double *stiffVecGen(double, double, int);
double randDouble(int, int);

int main (int argc, char *argv[]) {
    
    double pBond = 0.7, strain = 0.0, youngMod = 1.0;
    int netSize = 5;
    
    srand( time(NULL) );
    
    //Otherwise, parse the command line parameters.

    for (int i = 1; i < argc; i += 2) {
        string str = argv[i];
        if (!str.compare("-str")) {
            strain = atof(argv[i + 1]);
        } else if (!str.compare("-size")) {
            netSize = atoi(argv[i + 1]);
        } else if (!str.compare("-p")) {
            pBond = atof(argv[i + 1]);
        } else if (!str.compare("-y")) {
            youngMod = atof(argv[i + 1]);
        } else if (!str.compare("-help") || !str.compare("help")) {
            usageExit();
        } else {
            cout << argv[i] << " is an illegal argument. " << endl;
            usageExit();
        }
    }
    
    // Now that those are parsed, we can start to generate our network.
     
    Node ***nodeNetwork = new Node**[netSize];
    
    for (int i = 0; i < netSize; i++) {
        
        nodeNetwork[i] = new Node*[netSize];
        
        for (int j = 0; j < netSize; j++) {
            
            double *pos = new double [2];
            
            // x-coordinate
            pos[0] = RESTLEN * (i / 2.0  + j);
            
            // y-coordinate
            pos[1] = sqrt(3) / 2 * RESTLEN * i;
            
            double *sstiff = stiffVecGen(pBond, youngMod, 3);
            
            double rlen = RESTLEN;
            
            // Assign these data to the correct node
            
            nodeNetwork[i][j] = new Node(&*pos, &*sstiff, rlen);
        }
    }
    
    // TESTING CODE
    
    printf("Hurray! It all works!\n");
    printf("Now to test the energy function: \n");
    
    Funcd test;
    double energy = test(nodeNetwork, netSize, strain);
    printf("\nEnergy (drumroll, please): %5f\n", energy);
    
    

    /* AND THEN THIS......
    Frprmn<Funcd> frprmn(funcd);
    nodeNetwork = frprmn.minimalize(nodeNetwork);
    
    cout << funcd(nodeNetwork) << endl;
    */

	return 0;
}

//usageExit sends an error message about the usage of main and exits the 
//program.

void usageExit() {
    cout << "\nThe correct usage is " << endl;
    cout << "\n  main -str <strain> -size <network size> -p <bond ";
    cout << "probability> -y <young's modulus for spring>\n" << endl;
    exit(EXIT_FAILURE);
}

//stiffVecGen returns a vector containing three numbers, corresponding to the
//spring constants for the node. These numbers may be either yMod or 0.

double *stiffVecGen(double prob, double yMod, int numSprings) {
    
    double *ret = new double[numSprings];
    
    for (int i = 0; i < numSprings; i++)
        ret[i] = randDouble(0, 1) > prob ? 0 : yMod;

    return ret;
    
}

//randDouble returns a random number in the range [low, high).  
double randDouble(int low, int high) {
    
    double x = low + (high - low) * ((double) rand() / ((double) RAND_MAX + 1));
    return x;
    
}
