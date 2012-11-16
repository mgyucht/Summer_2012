/* frprmn.h
 * --------
 *
 * frprmn.h implements the Fletcher-Reeves-Polak-Ribiere method of minimizing
 * a multidimensional function using conjugate gradients.
 */

#ifndef FRPRMN_H_
#define FRPRMN_H_

#include "dlinmin.h"

extern int netSize;

template <class T>
struct Frprmn : Dlinemethod<T> {

    int iter;

    //Value of the function at the minimum.
    double fret;

    //Variables from a templated base class are not automatically inherited.

    using Dlinemethod<T>::func;
    using Dlinemethod<T>::linmin;
    using Dlinemethod<T>::p;
    using Dlinemethod<T>::xi;

    // Fractional tolerance.

    const double ftol;

    // Frprmn constructor.
    // Constructor arguments are funcd, the function or functor to be minimized, and an optional
    // argument ftoll, the fractional tolerance in the function value such that failure to decrease
    // by more than this amount on one iteration signals doneness.

    Frprmn(T &funcd, const double ftoll=1e-7) : Dlinemethod<T>(funcd),
    ftol(ftoll) {}

    // double *minimize(double *pp)
    // Given a starting point pp[0..n-1], performs the minimization on a function whose value
    // and gradient are provided by a functor funcd (see text).

    double *minimize(double *pp)
    {
        const int ITMAX = 1000000;
        const double EPS = 1.0e-18;
        const double GTOL = 1.0e-8;

        // Here ITMAX is the maximum allowed number of iterations; EPS is a small number to
        // rectify the special case of converging to exactly zero function value; and GTOL is the
        // convergence criterion for the zero gradient test.

        double gg, dgg;

        //Initializations.
        
        // n is the length of p.

        int n = 2 * netSize * netSize;
        
        p = pp;

        double *g = new double[n];
        double *h = new double[n];
        xi = new double[n];
        
        // fp is the energy of the network.
        // xi is the gradient of p.
        
        double fp = func(p);
        func.df(p, xi);
        
        // Sets g = h = -xi and negates xi.

        for (int j = 0; j < n; j++) {

            g[j] = -xi[j];
            xi[j] = h[j] = g[j];

        }

        // Loop over iterations.

        for (int its = 0; its < ITMAX; its++) {

            iter = its;
            
            fret = linmin(); //Next statement is one possible return:
            
            if (2.0 * abs(fret - fp) <= ftol * (abs(fret) + abs(fp) + EPS)) {
                
                delete[] g;
                delete[] h;
                delete[] xi;
                
                return p;

            }

            fp = fret;
            func.df(p, xi);

            //Test for convergence on zero gradient.
            double test = 0.0;
            double den = MAX(fp,1.0);

            for (int j = 0; j < n; j++) {

                double temp = abs(xi[j]) * MAX(abs(p[j]), 1.0) / den;

                if (temp > test)

                    test = temp;

            }

            // The other possible return.

            if (test < GTOL) {

                delete[] g;
                delete[] h;
                delete[] xi;
                
                #if DEBUGFRPRMN
                    
                    printf("Second return. Count is %d\n", its);
                    
                #endif
                
                return p;

            }

            dgg = gg = 0.0;

            for (int j = 0; j < n; j++) {

                gg += g[j]*g[j];

                // dgg += xi[j]*xi[j]; //This statement for Fletcher-Reeves.
                dgg += (xi[j]+g[j])*xi[j]; //This statement for Polak-Ribiere.

            }

            // Unlikely. If gradient is exactly zero, then we are already done.

            if (gg == 0.0) {

                delete[] g;
                delete[] h;
                delete[] xi;
                
                #if DEBUGFRPRMN
                    
                    printf("Third return. Count is %d\n", its);
                    
                #endif
                
                return p;

            }

            double gam = dgg/gg;

            for (int j = 0; j < n; j++) {

                g[j] = -xi[j];
                xi[j] = h[j] = g[j] + gam * h[j];

            }
        }
        throw("Too many iterations in frprmn");
    }
};

#endif /*FRPRMN_H_*/
