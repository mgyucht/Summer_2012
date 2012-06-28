#ifndef DBRENT_H_
#define DBRENT_H_

#include <limits>
#include "bracket.h"

// This struct implements Brent’s method to find a minimum, modified to use
// derivatives. Inherits the Bracketmethod class in bracket.h.

struct Dbrent : Bracketmethod {

    double xmin,fmin;
    const double tol;

    // Initializer for Dbrent. This sets the tolerance variable tol.

    Dbrent(const double toll=3.0e-8) : tol(toll) {}

    // Given a functor funcd that computes a function and also its derivative
    // function df, and given a bracketing triplet of abscissas ax, bx, cx
    // [such that bx is between ax and cx, and f(bx) is less than both f(ax)
    // and f(cx)], this routine isolates the minimum to a fractional precision
    // of about tol using a modification of Brent’s method that uses
    // derivatives. The abscissa of the minimum is returned as xmin, and the
    // minimum function value is returned as min, the returned function value.

    template <class T>
        double minimize(T &funcd) {

            const int ITMAX = 1000;
            const double ZEPS = std::numeric_limits<double>::epsilon()*1.0e-3;

            //Will be used as flags for whether pro-posed steps are acceptable or not.
            bool ok1,ok2;

            double a, b, d = 0.0, d1, d2, du, dv, dw, dx, e = 0.0;
            double fu, fv, fw, fx, olde, tol1, tol2, u, u1, u2, v, w, x, xm;

            // Comments following will point out only differences from the routine in Brent. Read
            // that routine first.

            // ax, bx, and cx are inherited from Bracketmethod. They define the
            // initial abcissa guesses.

            a = (ax < cx ? ax : cx);
            b = (ax > cx ? ax : cx);
            x = w = v = bx;
            fw = fv = fx = funcd(x);
            dw = dv = dx = funcd.df(x);

            // All our housekeeping chores are doubled by the necessity of moving
            // around derivative values as well as function values.

            for (double iter=0; iter < ITMAX; iter++) {

                xm = 0.5 * (a + b);
                tol2 = 2.0 * (tol1 = tol * abs(x) + ZEPS);

                // If the distance between x and xm is less than the minimum x
                // tolerance, we're done.

                if (abs(x - xm) <= (tol2 - 0.5 * (b - a))) {

                    fmin = fx;
                    return xmin = x;

                }

                if (abs(e) > tol1) {

                    //Initialize these d’s to an out-of-bracket value.

                    d1 = 2.0 * (b - a);
                    d2 = d1;

                    // Secant method with both points.
                    if (dw != dx)

                        d1=(w-x)*dx/(dx-dw);

                    if (dv != dx)

                        d2=(v-x)*dx/(dx-dv);

                    // Which of these two estimates of d shall we take? We will insist that they be
                    // within the bracket, and on the side pointed to by the derivative at x:

                    u1 = x + d1;
                    u2 = x + d2;

                    ok1 = (a - u1) * (u1 - b) > 0.0 && dx * d1 <= 0.0;
                    ok2 = (a - u2) * (u2 - b) > 0.0 && dx * d2 <= 0.0;

                    //Movement on the step before last.
                    olde = e;
                    e = d;

                    // Take only an acceptable d, and if both are acceptable, then take the smallest one.

                    if (ok1 || ok2) {

                        if (ok1 && ok2)

                            d = (abs(d1) < abs(d2) ? d1 : d2);

                        else if (ok1)

                            d = d1;

                        else

                            d = d2;

                        if (abs(d) <= abs(0.5 * olde)) {

                            u = x + d;

                            if (u - a < tol2 || b - u < tol2)

                                d = SIGN(tol1, xm - x);

                        } else {

                            // Decide which segment by the sign of the derivative.
                            d = 0.5 * (e = (dx >= 0.0 ? a - x : b - x));

                        }

                    } else {

                        d = 0.5 * (e = (dx >= 0.0 ? a - x : b - x));
                    }

                } else {

                    d = 0.5 * (e = (dx >= 0.0 ? a - x : b - x));

                }

                if (abs(d) >= tol1) {

                    u = x + d;
                    fu = funcd(u);

                } else {

                    u = x + SIGN(tol1, d);
                    fu = funcd(u);

                    // If the minimum step in the downhill direction takes us uphill, then we are done.

                    if (fu > fx) {

                        fmin=fx;
                        return xmin=x;

                    }
                }

                du=funcd.df(u);

                //Now all the housekeeping, sigh.

                if (fu <= fx) {

                    if (u >= x)

                        a=x;

                    else

                        b=x;

                    mov3(v,fv,dv,w,fw,dw);
                    mov3(w,fw,dw,x,fx,dx);
                    mov3(x,fx,dx,u,fu,du);

                } else {

                    if (u < x)

                        a=u;

                    else

                        b=u;

                    if (fu <= fw || w == x) {

                        mov3(v,fv,dv,w,fw,dw);
                        mov3(w,fw,dw,u,fu,du);

                    } else if (fu < fv || v == x || v == w) {

                        mov3(v,fv,dv,u,fu,du);

                    }

                }
            }
            throw("Too many iterations in routine dbrent");
        }
};

#endif /*DBRENT_H_*/
