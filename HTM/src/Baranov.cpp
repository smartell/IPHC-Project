// Baranov.cpp
#undef  MINSURV
#define MINSURV 0.01
double getFt(const double& m, const double& bt, double& ct);

dvector getFt(const double&m, const dvector& bt, dvector& ct)
{
	// Vectorized version
	int n = size_count(ct);
	double d_bt;
	double d_ct;
	double tmpf;
	dvector ft(1,n);
	for (int i = 1; i <= n; i++)
	{
		d_bt  = bt(i);
		d_ct  = ct(i);
		tmpf  = getFt(m,d_bt,d_ct);
		ft(i) = tmpf;
	}
	return(ft);
}

double getFt(const double& m, const double& bt, double& ct)
{
	/* |----------------------------------------------------------------------|
	   | Calculate instantaneous fishing mortality rate conditioned on catch
	   |----------------------------------------------------------------------|
	   | Args: M, natural mortality
	   |       bt, biomass at start of the year
	   |       ct, observed catch weight in the fishery.
	   |
	*/

	// initial guess for ft based on popes approximation
	double ctmp = ct;
	double ft,zt;

	ft = ctmp / (bt*exp(-0.5 * m));
	if(exp(-ft) < MINSURV)
	{
		ft   = -log(MINSURV);
		zt   = m + ft;
		ctmp = bt*(ft/zt*(1.0-exp(-zt)));

	}
	ct = ctmp;

	// iteratively solve the catch equaition
	double dct,ot,t0,t1,t2;
	for (int i = 0; i < 12; ++i)
	{	
		zt   = m + ft;
		t0   = exp(-zt);
		ot   = 1.0 - t0;
		ctmp = bt*(ft/zt*ot);

		// derivative of the catch equation
		t1   = ot/(zt*zt);
		t2   = t0/zt;
		dct  = bt*(-2.0*t1 + 2.0*t2 + 2.0*ft*t1/zt - 2.0*ft*t2/zt - ft*t2);

		// update Ft
		ft += (ctmp-ct)/dct;
		//cout<<i<<" "<<ft<<"\t"<<ctmp-ct<<endl;
	}
	return(ft);
}