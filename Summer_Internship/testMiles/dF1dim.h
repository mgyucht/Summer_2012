#ifndef DF1DIM_H_
#define DF1DIM_H_

template <class T>
struct Df1dim {
	const VecDoub &p;
	const VecDoub &xi;
	Int n;
	T &funcd;
	VecDoub xt;
	VecDoub dft;
	Df1dim(VecDoub_I &pp, VecDoub_I &xii, T &funcdd) : p(pp),
	xi(xii), n(pp.size()), funcd(funcdd), xt(n), dft(n) {}
	Doub operator()(const Doub x)
	{
		for (Int j=0;j<n;j++)
			xt[j]=p[j]+x*xi[j];
		return funcd(xt);
	}
	Doub df(const Doub x)
	{
		Doub df1=0.0;
		funcd.df(xt,dft);
		for (Int j=0;j<n;j++)
			df1 += dft[j]*xi[j];
		return df1;
	}
};


#endif /*F1dim_H_*/
