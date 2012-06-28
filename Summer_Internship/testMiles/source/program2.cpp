
#include <string>
#include <ctime>
#include <time.h>
#include <stdlib.h>
#include <math.h>
#include <fstream>
#include "funcd.h"
#include "frprmn.h"
#include "utils.h"
#include "debug.h"
#include "nonaffinity.h"

using namespace std;

int netSize = 20;
double strain;

int main() {
    
    double prob_bond = 0.8, strain = 0.01, youngs_mod = 1.0;
    
    double* positionArray = new double[2 * netSize * netSize];
    double*** springArray = new double**[netSize];
    
    for (int i = 0; i < netSize; i++) {
        springArray[i] = new double*[netSize];
        for (int j = 0; j < netSize; j++) {
            
            positionArray[(i * netSize + j) * 2] = j + i / 2.0 + i * strain 
                                                * sqrt(3.0) / 2.0;
            positionArray[(i * netSize + j) * 2 + 1] = sqrt(3.0) / 2.0 * i;
            
            springArray[i][j] = stiffVecGen(prob_bond, youngs_mod, 3);
        }
    }
    
    Funcd funcd(springArray);
    Frprmn<Funcd> frprmn(funcd);
    positionArray = frprmn.minimize(positionArray);
    double energy = funcd(positionArray);
    printf("%f\n", energy);

    return 0;
}
