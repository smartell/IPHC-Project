
/** \brief  A class object for a flexible growth model described by Schnute, 1981.
	
	Reference:
	Schnute, J. (1981). A versatile growth model with statistically stable 
	parameters. Canadian Journal of Fisheries and Aquatic Sciences, 
	38(9):1128–1140.
	
	
© Copyright 2012 International Pacific Halibut Commission - Martell. All Rights Reserved.

	\author Martell IPHC
	\author $LastChangedBy$
	\date 2012-10-20
	\date $LastChangedDate$
	\version 1	\sa
	
	The growth model is defined as follows:
	/f*
	l_j = ( y1^b + (y2^b-y1^b)*(1-exp(-a*(t-t1)))/(1-exp(-a*(t2-t1))) )^(1/b)
	*f/
	
	@param y1 length at t1
	@param y2 length at t2
	@param t1 reference time (age) for y1 length
	@param t2 reference time (age) for y2 length
	@param a  metabolic rate parameter (Brody growth coefficient)
	@param b  shape parameter (b=1 equates to von Bertalanffy growth)
	
**/

#ifndef _GROWTHMODEL_H_
#define _GROWTHMODEL_H_

class growthModel
{
private:
	double m_t1;
	double m_t2;
	double m_y1;
	double m_y2;
	double m_a;
	double m_b;
	
	
	dvector m_t;	//vector of ages (time)
	
public:
	// default constructor
	growthModel() { m_b = 1.0; }
	
	// destructor
	~growthModel() {}
	
	// member functions
	dvector length_at_age(	const dvector& t, 
							const double& y1, 
							const double& y2,
						  	const double& a, 
							const double& b,
							const double& t1,
							const double& t2);
	
};

#endif


dvector growthModel::length_at_age(	const dvector& t, 
									const double& y1, 
									const double& y2,
								  	const double& a, 
									const double& b,
									const double& j1,
									const double& j2)
{
	//l_j = ( y1^b + (y2^b-y1^b)*(1-exp(-a*(t-t1)))/(1-exp(-a*(t2-t1))) )^(1/b)
	// in log space
	cout<<y2<<endl;
	double  t1 = b*log(y1);
	double  t3 = exp(t1);
	double  t5 = exp(b*log(y2)) - t3;
	dvector t7 = (1.0-exp(-a*(t-j1))) / (1.0-exp(-a*(j2-j1)));
	double  t9 = 1.0/b;
	return ( exp(t9*log(t3 + t5 * t7)) );
	
}
