#ifndef DLINMIN_H_
#define DLINMIN_H_

#include "dbrent.h"
#include "dF1dim.h"


template <class T>
struct Dlinemethod {
//	Base class for line-minimization algorithms. Provides the line-minimization routine linmin.
	VecDoub p;
	VecDoub xi;
	T &func;
	Int n;
	Dlinemethod(T &funcc) : func(funcc) {}
//	Constructor argument is the user-supplied function or functor to be minimized.
	Doub linmin()
//	Line-minimization routine. Given an n-dimensional point p[0..n-1] and an n-dimensional
//	direction xi[0..n-1], moves and resets p to where the function or functor func(p) takes on
//	a minimum along the direction xi from p, and replaces xi by the actual vector displacement
//	that p was moved. Also returns the value of func at the returned location p. This is actually
//	all accomplished by calling the routines bracket and minimize of Brent.
	{
		Doub ax,xx,xmin;
		n=p.size();
		Df1dim<T> Df1dim(p,xi,func);
		ax=0.0; //Initial guess for brackets.
		xx=1.0;
		Dbrent Dbrent;
		Dbrent.bracket(ax,xx,Df1dim);
		xmin=Dbrent.minimize(Df1dim);
		for (Int j=0;j<n;j++) {// Construct the vector results to return.
			xi[j] *= xmin;
			p[j] += xi[j];
		}
		return Dbrent.fmin;
	}
	
	
};



#endif /*LINMIN_H_*/