#include <iostream>
#include <stdlib.h>
#include <time.h>
#include <math.h>

double randDouble(int low, int high);
double *stiffVecGen(double, double, int);

using namespace std;

int main () {
    
    srand( time(NULL) );
    
    cout << "Testing randDouble: \n" << endl;
    
    for (int i = 0; i < 5; i++)
        cout << randDouble(0, 1) << endl;
    
    cout << "Testing stiffVecGen: \n" << endl;
    
    for (int i = 0; i < 5; i++) {
        double * testArray = stiffVecGen(0.7, 0.8, 3);
        cout << "( ";
        for (int j = 0; j < 3; j++)
            cout << testArray[j] << " ";
        cout << ")" << endl;
    }
    
    return 0;
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
