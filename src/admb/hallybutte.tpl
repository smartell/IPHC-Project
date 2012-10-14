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
	// Model dimensions
	init_int narea;
	init_int ngroup;
	init_int nsex;
	init_int syr;
	init_int nyr;
	init_int sage;
	init_int nage;
	!! COUT(syr);
	!! COUT(nyr);
	!! COUT(nage);

PARAMETER_SECTION
	init_number log_bo;
	init_number log_steepness;
	init_number log_initR;
	init_number log_barR;
	
	init_vector log_m(1,nsex);
	
	objective_function_value f;
	
	number        bo;
	number steepness;
	
	vector     m(1,nsex);
	vector    pg(1,ngroup);
	
	5darray    N(1,narea,1,ngroup,1,nsex,syr,nyr,sage,nage);

INITIALIZATION_SECTION
	log_bo          7.0;
	log_steepness -0.30;
	log_m         -1.90;


PROCEDURE_SECTION
	/*    MAIN FUNCTION CALLS    */
	
	initParameters();
	
	initPopulationModel();
	
	/* ------------------------- */

	cout<<"Hello World"<<endl;
	N.initialize();
	exit(1);


FUNCTION initParameters
  {
	/** Transform parameters */
	bo        = mfexp(log_bo);
	steepness = mfexp(log_steepness);
	m         = mfexp(log_m);
	
	/* normal distribution for group proportions */
	dvector x(1,ngroup);
	x.fill_seqadd(1,1);
	x  = 2.*(-1. + 2.*(x-1.)/(ngroup-1.));
	pg = 1.0/sqrt(2.0*PI)*exp(-0.5*square(x));
	

	
	
  }

FUNCTION initPopulationModel
  {
	/** Initiaize numbers at age */
	
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

