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
#include <time.h>
#include "nr3.h"
#include "funcd.h"
#include <math.h>
#include "node.h"

using namespace std;

vector<double> stiffVecGen(double, double, int);
double randDouble(int, int);

int main (int argc, char *argv[]) {
    
    double strain = 0.0, pBond = 1.0, youngMod = 1.0;
    int netSize = 200;
    srand( time(NULL) );
    
    // If there are not the right number of arguments, end the program's
    // execution.

    if (argc != 9) {
        cout << "The correct usage is " << endl;
        cout << "\n main -str <strain> -size <network size> -p <bond ";
        cout << "probability> -y <young's modulus for spring>" << endl;
        exit(0);
    }
    
    //Otherwise, parse the command line parameters.

    for (int i = 1; i < 9; i += 2) {
        if (argv[i] == "-str")
            strain = atof(argv[i + 1]);
        else if (argv[i] == "-size")
            netSize = atoi(argv[i + 1]);
        else if (argv[i] == "p")
            pBond = atof(argv[i + 1]);
        else if (argv[i] == "y")
            youngMod = atof(argv[i + 1]);
        else {
            cout << argv[i] << " is an illegal argument. " << endl;
            exit(0);
        }
    }
    
    // Now that those are parsed, we can start to generate our network.
     
    VecNode2D nodeNetwork(netSize, netSize);

    for (int i = 0; i < netSize; i++) {
        for (int j = 0; j < netSize; j++) {
            
            double pos[] = { i * RESTLEN, j * RESTLEN };
            vector<double> position (pos, pos + sizeof(pos) / sizeof(double));
            vector<double> initLen (3, RESTLEN);
            vector<double> sstiff = stiffVecGen(pBond, youngMod, 3);
            
            // Assign this vector to the correct node
        }
    }
    
    // TESTING CODE
    
    for (int i = 0; i < netSize; i++) {
        for (int j = 0; j < netSize; j++) {
            cout << " " << nodeNetwork[i][j].position[1] << \
                " " << nodeNetwork[i][j].position[2];
        }
        cout << endl;
    }

    /* AND THEN THIS......
    Funcd funcd;
    Frprmn<Funcd> frprmn(funcd);
    nodeNetwork = frprmn.minimalize(nodeNetwork);
    
    cout << funcd(nodeNetwork) << endl;
    */
    
	return 0;
}

//stiffVecGen returns a vector containing three numbers, corresponding to the
//spring constants for the node. These numbers may be either yMod or 0.

vector<double> stiffVecGen(double prob, double yMod, int numSprings) {
    double array[numSprings];
    
    for (int i = 0; i < numSprings; i++) {
        array[i] = randDouble(0, 1) > prob ? 0 : yMod;
    }

    vector<double> result (array, array + sizeof(array) / sizeof(double));
    return result;
}

//randDouble returns a random number in the range [low, high). 

double randDouble(int low, int high) {
    double x = low + (high - low) * ((double) rand() / ((double) RAND_MAX + 1));
    return x;
}
