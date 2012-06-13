// dlinmin.h
// ---------
//
// dlinmin.h contains the struct Dlinemethod, which uses a variant of Brent's
// method to perform a line minimization for a multidimensional function.

#ifndef DLINMIN_H_
#define DLINMIN_H_

#include "dbrent.h"
#include "dF1dim.h"

// Base class for line-minimization algorithms. Provides the line-minimization
// routine linmin.

template <class T>
struct Dlinemethod {
    
	VecDoub p;
	VecDoub xi;
	T &func;
	Int n;
    
    // Constructor argument is the user-supplied function or functor to be 
    // minimized.
    
	Dlinemethod(T &funcc) : func(funcc) {}
    
    // Line-minimization routine. Given an n-dimensional point p[0..n-1] and
    // an n-dimensional direction xi[0..n-1], moves and resets p to where the
    // function or functor func(p) takes on a minimum along the direction xi from
    // p, and replaces xi by the actual vector displacement that p was moved. Also
    // returns the value of func at the returned location p. This is actually
    // all accomplished by calling the routines bracket and minimize of Brent.
    
	Doub linmin() {
        
		Doub ax, xx, xmin;
		n = p.size();
		Df1dim<T> Df1dim(p,xi,func);
        
        //Initial guess for brackets.
        
		ax=0.0; 
		xx=1.0;
        
		Dbrent Dbrent;
		Dbrent.bracket(ax,xx,Df1dim);
		xmin=Dbrent.minimize(Df1dim);
        
        // Construct the vector results to return.
        
		for (Int j=0;j<n;j++) {
            
			xi[j] *= xmin;
			p[j] += xi[j];
            
		}
        
		return Dbrent.fmin;
	}
};



#endif /*LINMIN_H_*/
