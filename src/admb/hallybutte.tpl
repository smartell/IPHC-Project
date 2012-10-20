//  ***********************************************************************  //
//  hallybutte (from middle english meaning "holy and flat")                 //
//                                                                           //
//  Created by Martell on 2012-10-12.                                        //
//  Copyright (c) 2012. All rights reserved.                                 //
//  Comments:                                                                //
//                                                                           //
//                                                                           //
//                                                                           //
//                                                                           //
//                                                                           //
//                                                                           //
//  Indexs:                                                                  //
//         area   f                                                          //
//         group  g                                                          //
//         sex    h                                                          //
//         year   i                                                          //
//         age    j                                                          //
//         grp    k                                                          //
//                                                                           //
//  Array dimensions:                                                        //
//         collapse area group and sex into a single dimension, then create  //
//         a linked-list to index sex area and group.                        //
//                                                                           //
//                                                                           //
//  MOVEMENT MODEL BASED ON GRAVITY WEIGTHS AND RESIDENCY                    //
//  -reference is Carruthers et al. 2011 and Caddy 1975                      //
//  -estimate a vector of area-specific gravity weights and a residencey     //
//   parameter (probability of staying put in all areas).                    //
//  - set M_{f,f'} = r+g_{f}  if f = f'                                      //
//    set M_{f,f'} = g_{f}    if f!= f'                                      //
//        M_{f,f'} = M_{f,f'}/sum_{f'} M_{f,f'}                              //
//  - M_{f,f'} is the probability of movement from f to f'                   //
//                                                                           //
//                                                                           //
//                                                                           //
//  ***********************************************************************  //


DATA_SECTION
	// Counters
	int f;		// area  
	int g;      // group 
	int h;      // sex   
	int i;      // year  
	int j;      // age   
	int k;      // grp   
	
	// Model dimensions
	init_int narea;    
	init_int ngroup;   
	init_int nsex;     
	init_int syr;      
	init_int nyr;      
	init_int sage;     
	init_int nage;     
	!! ECHO(narea )
	!! ECHO(ngroup)
	!! ECHO(nsex  )
	!! ECHO(syr   )
	!! ECHO(nyr   )
	!! ECHO(sage  )
	!! ECHO(nage  )
	ivector i_age(sage,nage);
	vector  d_age(sage,nage);
	!! i_age.fill_seqadd(sage,1);
	!! d_age.fill_seqadd(sage,1);
	
	// Linked list to manage arrays
	int N_grp;
	!!  N_grp = narea*ngroup*nsex;
	ivector   i_sex(1,N_grp);
	ivector  i_area(1,N_grp);
	ivector i_group(1,N_grp);
	LOC_CALCS
		k = 0;
		for(f=1;f<=narea;f++)
			for(g=1;g<=ngroup;g++)
				for(h=1;h<=nsex;h++)
				{
					k ++;
					i_sex(k)   = h;
					i_area(k)  = f;
					i_group(k) = g;
				}	
		COUT(i_sex);
		COUT(i_area);
		COUT(i_group);
	END_CALCS
	
	
	// Growth and allometry parameters
	init_vector l_sage(1,nsex);
	init_vector l_nage(1,nsex);
	init_vector    vbk(1,nsex);
	init_vector      a(1,nsex);
	init_vector      b(1,nsex);
	!! ECHO(l_sage)
	!! ECHO(l_nage)
	!! ECHO(   vbk)
	!! ECHO(     a)
	!! ECHO(     b)
	
	// Female maturity at age
	init_number mat_age50;
	init_number mat_std50;
	!! ECHO(mat_age50)
	!! ECHO(mat_std50)
	

PARAMETER_SECTION
	init_number log_bo;
	init_number log_steepness;
	init_number log_initR;
	init_number log_barR;
	
	init_vector log_m(1,nsex);
	
	init_bounded_dev_vector log_initR_devs(sage,nage,-15.0,15.0,2);
	init_bounded_dev_vector log_barR_devs(syr+1,nyr,-15.0,15.0,2);
	
	objective_function_value f;
	
	number        bo;
	number steepness;
	
	vector     m(1,nsex);
	vector    pg(1,ngroup);
	
	3darray    N(syr,nyr,1,N_grp,sage,nage);

INITIALIZATION_SECTION
	log_bo          7.0;
	log_steepness -0.30;
	log_m         -1.90;


PROCEDURE_SECTION
	/*    MAIN FUNCTION CALLS    */
	
	initParameters();
	
	calcGrowth();
	
	initPopulationModel();
	
	runPopulationModel();
	
	/* ------------------------- */
	exit(1);


FUNCTION initParameters
  {
	/** Transform parameters */
	bo        = mfexp(log_bo);
	steepness = mfexp(log_steepness);
	m         = mfexp(log_m);
	
	/* normal distribution for group proportions */
	if( ngroup >1 )
	{
		dvector x(1,ngroup);
		x.fill_seqadd(1,1);
		x  = 2.*(-1. + 2.*(x-1.)/(ngroup-1.));
		pg = 1.0/sqrt(2.0*PI)*exp(-0.5*square(x));
		pg/= sum(pg);		
	}
	else
	{
		pg = 1.0;
	}
	
  }

FUNCTION calcGrowth
  {
	/** Calculate length-at-age and weight-at-age */
	growthModel c_growth;
	cout<<"GROWTH"<<endl;
	cout<<c_growth.length_at_age(d_age,10,100,0.2,1,sage,20)<<endl;
  }

FUNCTION initPopulationModel
  {
	
	N.initialize();
	
	dvariable tmpR;
	for(k=1;k<=N_grp;k++)
	{
		f = i_area(k);
		g = i_group(k);
		h = i_sex(k);
		
		/** Initiaize numbers at age in syr */
		for(j=sage;j<=nage;j++)
		{
			tmpR         = mfexp( log_initR + log_initR_devs(j) );
			tmpR        /= (nsex * narea); //dispersal kernel here.
			N(syr)(k)(j) = pg(g) * tmpR * mfexp(-m(h)*(j-sage));
			if( j==nage )
			{
				N(syr)(k)(j) /= ( 1.0 - mfexp(-m(h)) );
			}
		}
		
		/** Initial recruitment in each year */
		for(i=syr+1;i<=nyr;i++)
		{
			tmpR          = mfexp( log_barR + log_barR_devs(i) );
			tmpR         /= (nsex * narea);	//dispersal kernel here.
			N(i)(k)(sage) = pg(g) * tmpR;
		}	
	}
	
  }

FUNCTION runPopulationModel
  {
	/** Update numbers-at-age in all N_grp and years */
	
	
	for(k=1;k<=N_grp;k++)
	{
		f = i_area(k);
		g = i_group(k);
		h = i_sex(k);
		for(i=syr+1;i<=nyr;i++)
		{
			N(i)(k)(sage+1,nage) =++ N(i-1)(k)(sage,nage-1)*exp(-m(h));
			N(i)(k)(nage)       +=   N(i-1)(k)(nage)*exp(-m(h));
		}
	}
	//COUT(N)
  }

REPORT_SECTION


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
	#define COUT(object) cout << #object "\n" << object <<endl;
	
	#undef ECHO
	#define ECHO(object) echoinput << "# " #object "\n" << object <<endl;

	#include <admodel.h>
	#include <time.h>
	#include <statsLib.h>
	#include <growth.cpp>
	time_t start,finish;
	long hour,minute,second;
	double elapsed_time;
	
	ofstream echoinput("echoinput.txt");
	
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

