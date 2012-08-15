// dF1dim.h
// --------
//
// dF1dim.h provides the struct used in the linmin() method provided in
// class Dlinemethod in dlinmin.h.

#ifndef DF1DIM_H_
#define DF1DIM_H_

extern int netSize;

template <class T>
struct Df1dim {

    // p is the position array. xi is the gradient array.
    
    const double *p;
    const double *xi;
    T &funcd;
    double *xt;
    double *dft;

    Df1dim(double *pp, double *xii, T &funcdd) :
        p(pp), xi(xii), funcd(funcdd) {

        xt = new double[2 * netSize * netSize];
        dft = new double[2 * netSize * netSize];
        
    }
    
    ~Df1dim() {
        
        delete[] xt;
        delete[] dft;
    
    }
        
    double operator() (const double x) {

        for (int j = 0; j < 2 * netSize * netSize; j++)

            xt[j]=p[j]+x*xi[j];

        return funcd(xt);

    }

    double df(const double x) {

        double df1 = 0.0;
        funcd.df(xt, dft);

        for (int j = 0; j < 2 * netSize * netSize; j++)

            df1 += dft[j] * xi[j];

        return df1;
    }
};

#endif /*DF1DIM_H_*/
