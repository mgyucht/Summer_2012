/* batch.cpp
 * --------
 *
 * Author: Miles Yucht
 * Date: Mon June 27 2012
 */

#include <string>
#include <stdlib.h>
#include <stdio.h>

using namespace std;

int main (int argc, char *argv[]) {
    
    for (double pBond = 0.9; pBond > 0.5; pBond -= 0.025) {
        FILE * file = fopen("compare_data.txt", "a");
        fprintf(file, "pBond = %4.3f\n", pBond); 
        fclose(file);
        for (double strain = 0.01; strain < 0.095; strain += 0.01) {
            
            char command[32];
            sprintf(command, "./program -str %3.2f -p %4.3f", strain, pBond);
            system(command);

        }
    }

    return 0;
}
