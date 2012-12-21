//  ******************************************************************
//  | HTM Halibut Toy Model                                          |
//  |                                                                |
//  | Created by Martell on 2012-11-16.                              |
//  | Copyright (c) 2012. All rights reserved.                       |
//  | Comments: A spatial delay difference model                     |
//  | See the ReadMe.md file for more information about this program.|
//  |                                                                |
//  ******************************************************************


DATA_SECTION
	
	int simFlag;
	int rseed;
	int retroYears;
	LOC_CALCS
		simFlag=0;
		rseed=999;
		int on,opt;
		//the following line checks for the "-SimFlag" command line option
		//if it exists the if statement retreives the random number seed
		//that is required for the simulation model
		if((on=option_match(ad_comm::argc,ad_comm::argv,"-sim",opt))>-1)
		{
			simFlag = 1;
			rseed   = atoi(ad_comm::argv[on+1]);
			//if(SimFlag)exit(1);
		}
		
		// command line option for retrospective analysis. "-retro retro_yrs"
		retroYears=0;
		if((on=option_match(ad_comm::argc,ad_comm::argv,"-retro",opt))>-1)
		{
			retroYears=atoi(ad_comm::argv[on+1]);
			cout<<"______________________________________________________\n"<<endl;
			cout<<"    **Implementing Retrospective analysis** "<<endl;
			cout<<"    **Number of retrospective years = "<<retroYears<<endl;
			cout<<"______________________________________________________"<<endl;
		}
	END_CALCS
	
	
	// Model dimensions
	init_int syr;
	init_int nyr;
	init_int agek;
	init_int narea;
	
	// Growth parameters
	init_number alpha;
	init_number rho;

	// Total removals by year and area
	init_matrix removals(syr,nyr,0,narea);
	matrix 		  obs_ct(syr,nyr,1,narea);
	LOC_CALCS
		obs_ct = trans(trans(removals).sub(1,narea));
	END_CALCS

	// Setline survey WPUE
	init_int	n_surv_obs;
	init_matrix survey_wpue(1,n_surv_obs,0,narea);
	
	// Commercial fishery WPUE
	init_int 	n_fish_obs
	init_matrix fishery_wpue(1,n_fish_obs,0,narea);

	// Movement matrix
	init_matrix 		moveProb(1,narea,1,narea);
	vector 	recruitmentDispersal(1,narea);
	LOC_CALCS
		int i;
		/* Ensure each row sums to 1*/
		for(i=1;i<=narea;i++)
		{
			moveProb(i) = moveProb(i)/sum(moveProb(i));
		}
	END_CALCS


INITIALIZATION_SECTION
	log_bo    7.0;
	steepness 0.7;
	m         0.15;
	log_rbar  1.122648;

PARAMETER_SECTION
	init_number log_bo;
	init_bounded_number steepness(0.2,1.0);
	init_bounded_number m(0.01,0.4,-2);
	init_number log_rbar;
	init_bounded_dev_vector log_rbar_devs(syr,nyr,-5,5,2);
	
	init_vector log_fbar(1,narea);
	init_bounded_matrix log_f_devs(syr,nyr,1,narea,-5,5,2);
	
	init_vector log_q(1,narea);
	
	objective_function_value f;
	
	number bo;
	number ro;
	number no;
	number reck;
	number rbar;
	number wbar;
	number wk;
	number so;
	number beta;
	number fpen;
	
	vector p(1,narea);
	vector q(1,narea);
	
	matrix rt(syr,nyr,1,narea);
	matrix ft(syr,nyr,1,narea);
	matrix ct(syr,nyr,1,narea);
	matrix wt(syr,nyr,1,narea);
	matrix it(syr,nyr,1,narea);
	matrix yt(syr,nyr,1,narea);
	
	matrix bt(syr,nyr+1,1,narea);
	matrix nt(syr,nyr+1,1,narea);
	
	

PRELIMINARY_CALCS_SECTION
	/* calculate dispersal vector */
	initializeDispersalKernel(recruitmentDispersal,moveProb);
	
	/* if simFlag is true then generate simulated data */
	if(simFlag)
	{
		cout<<"|______________________________________________________|"<<endl;
		cout<<"| Simulation: overwriting data with simulated data     |\n";
		cout<<"|______________________________________________________|"<<endl;
		cout<<"| Not fully implemented yet.."<<endl;
		exit(1);
	}


PROCEDURE_SECTION
	initialize_model();
	
	calcFishingMortality();
	
	population_dynamics();
	
	calcFisheryCatchStatistics();
	
	calcSurveyStatistics();
	
	
FUNCTION initialize_model
	int i;
	dvariable s;
	
	s    = mfexp(-m);
	bo   = mfexp(log_bo);
	reck = 4.0* steepness / (1.0 - steepness);
	rbar = mfexp(log_rbar);
	wk   = alpha*( 1.0 - pow(rho,agek) )/(1.0-rho);
	wbar = (s*alpha + wk*(1.0-s))/(1.0-rho*s);
	no   = bo/wbar;
	ro   = no*(1.0-s);
	
	so   = reck*ro/bo;
	beta = (reck-1.0)/bo;
	
	
	/*
	|----------------------------------------------------------------------------
	| Initialize recruitment dispersal kernel (p)
	| - p is a vector of proportions of total recruitment that end up in a given 
	|   regulatory area.  
	| 
	| - This vector is based on the eigen vector of the movement probability with 
	|   an eigen value = 1.0.
	| 
	| - Could also add a vector of estimated deviates to represent the notion
	|   that dispersal of fish less than age-k is incomplete after age-k.
	|----------------------------------------------------------------------------
	*/
	p           = recruitmentDispersal;
	// COUT(recruitmentDispersal);

	// Initialize biomass and numbers;
	bt(syr)     = bo * p;
	nt(syr)     = no * p;
	for(i=syr;i<=nyr;i++)
	{
		rt(i)   = mfexp(log_rbar + log_rbar_devs(i)) * p;
		rt(i)   = ro*p;
	}
	// COUT(rt);
	
FUNCTION calcFishingMortality
  {
	/*
		| This routine fills the matrix of fishing mortality rates
		| in each area and ensues each column of f_devs sums to zero.
	*/
	
	int i,j;
	
	for(i=syr;i<=nyr;i++)
	{
		ft(i) = 1.0 - mfexp( log_fbar + log_f_devs(i) );
	}
	
	// ft.sub(2000,2009) = 0.4;
	// ft(1996,4) = 0.5;
	
	/* Penalty for log_f_devs */
	dvariable s;
	fpen.initialize();
	for(j=1;j<=narea;j++)
	{
		s     = mean(column(log_f_devs,j));
		fpen += 1000. * s*s;
	}
  }
	
FUNCTION void initializeDispersalKernel(dvector& d, const dmatrix& M)
  {
	/*
		|Call from PRELIMINARY_CALCS_SECTION only.
		|
		| This routine initializes the dispersal kernel for recruitment.
		| Numerically compute the dispersal vector recruitment_dispersal
		| by successive multiplication of the movement matrix (M)
		|
		| This method is pretty inefficient and I would not use this for
		| dvariable calculations.  
	*/
	int iter;
	int imin = d.indexmin();
	int imax = d.indexmax();
	dvector p(imin,imax);
	dvector pminus(imin,imax);
	p = 1.0/(imax-imin+1);
	
	for(iter=1;iter<=2000;iter++)
	{
		pminus = p;
		p = p * M;
		
		// cout<<iter<<" "<<norm2(p-pminus)<<endl;
		if(norm2(p-pminus)<=1.e-20) break;
	}
	d = p;
	return;
  }


FUNCTION  population_dynamics
  {
	int i,j;
	dvar_vector sj(1,narea);
	for(i=syr;i<=nyr;i++)
	{
		sj      = mfexp( -m - ft(i) );
		nt(i+1) = ( elem_prod( sj, nt(i) ) + rt(i)) * moveProb;
		bt(i+1) = ( elem_prod( sj, alpha*nt(i)+rho*bt(i) ) + wk*rt(i) ) * moveProb;
	}
	// COUT(bt);
  }
	
FUNCTION calcFisheryCatchStatistics
  {
	/*
		| Compute the predicted fishery catch statistics:
		| 1) commercial catch by area          (ct)
		| 2) commercial wpue by area           (it)
		| 3) mean weight of the catch by area  (wt)
	*/
	
	int i;
	dvar_vector z(1,narea);
	
	q = mfexp(log_q);
	for(i=syr;i<=nyr;i++)
	{
		/*1) commercial catch by area */
		z     = m + ft(i);
		ct(i) = elem_prod(elem_prod(elem_div(ft(i),z),1.-exp(-z)),bt(i));
		
		/*2) commercial wpue */
		it(i) = elem_prod(q,bt(i));
		
		/*3) mean weight of the catch */
		wt(i) = elem_div(bt(i),nt(i));
	}
	
  }

FUNCTION calcSurveyStatistics
  {
	/*
		| Compute predicted survey observations:
		| 1) relative abundance index for survey (yt)
	*/
	
	int i;
	dvar_vector qs(1,narea);
	qs = 1.0;
	
	for(i=syr;i<=nyr;i++)
	{
		/*1) relative abundance index for survey (yt) */
		yt(i)  = qs * bt(i);
	}
	
  }
	
REPORT_SECTION
	REPORT(nt);
	REPORT(bt);


TOP_OF_MAIN_SECTION
	time(&start);
	arrmblsize = 50000000;
	gradient_structure::set_GRADSTACK_BUFFER_SIZE(1.e7);
	gradient_structure::set_CMPDIF_BUFFER_SIZE(1.e7);
	gradient_structure::set_MAX_NVAR_OFFSET(5000);
	gradient_structure::set_NUM_DEPENDENT_VARIABLES(5000);
 

GLOBALS_SECTION
	/**
	\def REPORT(object)
	Prints name and value of \a object on ADMB report %ofstream file.
	*/
	#undef REPORT
	#define REPORT(object) report << #object "\n" << object << endl;
	
	#undef COUT
	#define COUT(object) cout<<fixed<<#object "\n"<<object<<endl;
	
	#include <iostream>
	#include <iomanip>
	using namespace std;
	
	#include <admodel.h>
	#include <time.h>
	#include <statsLib.h>
	time_t start,finish;
	long hour,minute,second;
	double elapsed_time;
	
FINAL_SECTION
	time(&finish);
	elapsed_time=difftime(finish,start);
	hour=long(elapsed_time)/3600;
	minute=long(elapsed_time)%3600/60;
	second=(long(elapsed_time)%3600)%60;
	cout<<endl<<endl<<"*******************************************"<<endl;
	cout<<"--Start time: "<<ctime(&start)<<endl;
	cout<<"--Finish time: "<<ctime(&finish)<<endl;
	cout<<"--Runtime: ";
	cout<<hour<<" hours, "<<minute<<" minutes, "<<second<<" seconds"<<endl;
	cout<<"*******************************************"<<endl;

