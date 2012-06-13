#ifndef BRACKET_H_
#define BRACKET_H_

//	Base class for one-dimensional minimization routines. Provides a routine to bracket a minimum
//	and several utility functions.

struct Bracketmethod { 
    
	Doub ax,bx,cx,fa,fb,fc;
    
	// Given a function or functor func, and given distinct initial
	// points ax and bx, this routine searches in the downhill direction
    // (defined by the function as evaluated at the initial points) 
    // and returns new points ax, bx, cx that bracket a minimum of the 
    // function. Also returned are the function values at the three 
    // points, fa, fb, and fc.
    
	template <class T>
	void bracket(const Doub a, const Doub b, T &func)
	{
        //	Here GOLD is the default ratio by which successive intervals are magnified and GLIMIT
        //	is the maximum magnification allowed for a parabolic-fit step.
        
		const Doub GOLD=1.618034,GLIMIT=100.0,TINY=1.0e-20;
        
		ax = a; 
        bx = b;
		Doub fu;
		fa = func(ax);
		fb = func(bx);
        
        //Switch roles of a and b so that we can go
		// downhill in the direction from a to b.
        
        if (fb > fa) { 
			SWAP(ax,bx);
			SWAP(fb,fa);
		}
        
        //First guess for c.
        
		cx = bx + GOLD * (bx - ax); 
        fc = func(cx);
        
        //Keep returning here until we bracket.
        
        while (fb > fc) { 
            
            // Compute u by parabolic extrapolation from a, b, c. TINY is used 
            // to prevent any possible division by zero.
            
            Doub r = (bx - ax) * (fb - fc); 
            Doub q = (bx - cx) * (fb - fa);
            Doub u = bx - ((bx - cx) * q - (bx - ax) * r) / (2.0 * SIGN(MAX(abs(q - r), TINY), q - r));
            Doub ulim = bx + GLIMIT * (cx - bx);
            
            //We won’t go farther than this. Test various possibilities:
            
            //Parabolic u is between b and c: try it.
            if ((bx-u)*(u-cx) > 0.0) { 
                
                fu=func(u);
                
                if (fu < fc) {
                    // Got a minimum between b and c.
                    
                    ax=bx;
                    bx=u;
                    fa=fb;
                    fb=fu;
                    return;
                } else if (fu > fb) { 
                    //Got a minimum between between a and u.
                    
                    cx=u;
                    fc=fu;
                    return;
                }

                //Parabolic fit was no use. Use default magfnification.

                u=cx+GOLD*(cx-bx); 
                fu=func(u); 

            } else if ((cx-u)*(u-ulim) > 0.0) { 
                //Parabolic fit is between c and its allowed limit.

                fu=func(u); 
                if (fu < fc) {
                    shft3(bx,cx,u,u+GOLD*(u-cx));
                    shft3(fb,fc,fu,func(u));
                }
            } else if ((u-ulim)*(ulim-cx) >= 0.0) { 
                //Limit parabolic u to maximum allowed value.

                u=ulim; 
                fu=func(u);
            } else {
                // Reject parabolic u, use default magnificaution.

                u=cx+GOLD*(cx-bx); 
                fu=func(u);
            }

            // Eliminate oldest point and continue.
            
            shft3(ax,bx,cx,u);
            shft3(fa,fb,fc,fu);
        }
    }
    
    //	Utility function used in this structure or others derived from it.
    inline void shft2(Doub &a, Doub &b, const Doub c)
	{
		a=b;
		b=c;
	}
    
	inline void shft3(Doub &a, Doub &b, Doub &c, const Doub d)
	{
		a=b;
		b=c;
		c=d;
	}
    
	inline void mov3(Doub &a, Doub &b, Doub &c, const Doub d, const Doub e,
					 const Doub f)
	{
		a=d; 
        b=e; 
        c=f;
	}
};

#endif /*BRACKET_H_*/
