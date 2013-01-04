//  ******************************************************************
//  | HTM Halibut Toy Model                                          |
//  |                                                                |
//  | Created by Martell on 2012-11-16.                              |
//  | Copyright (c) 2012. All rights reserved.                       |
//  | Comments: A spatial delay difference model                     |
//  | See the ReadMe.md file for more information about this program.|
//  |                                                                |
//  |                                                                |
//  |                                                                |
//  |                                                                |
//  |  TODO: 
//  |      - Add stock recruitment function  DONE 
//  |      - Add recruitemnt residuals to objective function. DONE
//  |      - Random walk in growth. DONE
//  |      - Simulation model, based on S_R relationship.
//  |                                                                |
//  |                                                                |
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
	init_number d_alpha;
	init_number rho;

	// Habitat area
	init_vector habitatArea(1,narea);
	!! habitatArea = habitatArea / 1.e5;

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

	// Average weight of the commercial catch
	init_int n_wbar_obs
	init_matrix fishery_wbar(1,n_wbar_obs,0,narea);

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
	log_fbar -2.302585;

PARAMETER_SECTION
	init_number log_bo(1);
	init_bounded_number steepness(0.2,1.0);
	init_bounded_number m(0.01,0.4,-2);
	init_bounded_vector log_rbar(1,narea,-9,9,1);
	init_bounded_matrix log_rbar_devs(syr,nyr,1,narea,-15,15,2);
	
	// fishing mortality
	init_vector log_fbar(1,narea);
	init_bounded_matrix log_f_devs(syr,nyr,1,narea,-5,5,2);
	
	// gravity model parameters
	init_bounded_number log_residency(-5,5,3);
	init_bounded_dev_vector log_gravity_weights(1,narea,-5,5,3);

	// catchability parameters for commercial fishery and survey
	init_number  log_q;
	init_bounded_dev_vector log_q_dev(1,narea,-5,5,2);
	init_number log_qs;
	init_bounded_dev_vector log_qs_dev(1,narea,-5,5,2);

	// random walk parameters for growth intercept.
	init_number log_alpha(1);
	init_bounded_vector log_alpha_dev(syr,nyr-1,-5,5,4);
	
	objective_function_value f;
	
	number bo;
	number ro;
	number no;
	number reck;
	number wbar;
	number wk;
	number so;
	number beta;
	number Delta;							// Residencey parameter
	number fpen;
	number rpen;
	
	vector            p(1,narea);
	vector            q(1,narea);
	vector           qs(1,narea);
	vector         rbar(1,narea);
	vector      gravity(1,narea);			// Gravity weights
	vector    ell_catch(1,narea);
	vector ell_wpue_fsh(1,narea);
	vector ell_wbar_fsh(1,narea);
	vector ell_wpue_srv(1,narea);
	vector       ell_Pj(1,narea);

	vector    St(syr,nyr);					// Spawning biomass
	vector    Rt(syr,nyr);					// Recruitment numbers
	vector alpha(syr,nyr);					// Growth intercept vector.
	vector  Rhat(syr+agek,nyr+agek);
	vector delta(syr+agek,nyr);
	
	matrix     rt(syr,nyr,1,narea);
	matrix     ft(syr,nyr,1,narea);
	matrix     ct(syr,nyr,1,narea);
	matrix     wt(syr,nyr,1,narea);
	matrix     it(syr,nyr,1,narea);
	matrix     yt(syr,nyr,1,narea);
	matrix varphi(syr,nyr,1,narea);			// Catch residuals

	matrix 	   eta(1,n_fish_obs,1,narea);	// Commercial WPUE residuals
	matrix  	nu(1,n_wbar_obs,1,narea);	// Commercial wbar residuals
	matrix epsilon(1,n_surv_obs,1,narea);	// Setline survey residuals WPUE
	
	matrix bt(syr,nyr+1,1,narea);
	matrix nt(syr,nyr+1,1,narea);

	matrix Pj(1,narea,1,narea);
	
	

PRELIMINARY_CALCS_SECTION
	/* calculate dispersal vector */
	//initializeDispersalKernel(recruitmentDispersal,moveProb);
	
	/* if simFlag is true then generate simulated data */
	if(simFlag)
	{
		cout<<"|______________________________________________________|"<<endl;
		cout<<"| Simulation: overwriting data with simulated data     |\n";
		cout<<"|______________________________________________________|"<<endl;
		cout<<"| Not fully implemented yet.."<<rseed<<endl;

		runSimulationModel(rseed);
		cout<<"|"<<endl;
		
	}


PROCEDURE_SECTION

	calcFishingMortality();

	growthModel();
	
	initializeModel();

	calcGravityModel();
	
	initialStates();

	population_dynamics();
	
	calcFisheryCatchStatistics();
	
	calcSurveyStatistics();

	calcStockRecruit();

	calcObjectiveFunction();

FUNCTION void runSimulationModel(int& seed)
  {
	/* 
	| -------------------------------------------------------------------- |
	| Simulation model conditioned on the input catch data and dispersel
	| kernel based on tagging data.  The purpose of this simulation model
	| is to determine how estimable the model parameters are, given WPUE
	| data only.
	| -------------------------------------------------------------------- |
	| 
	| Simulation model is based on a stock-recruitment relationship. 
	| Must calculate total recruitment and apportion it among regulatory 
	| areas, such that average fishing mortality rates are achieved based
	| based on recruitment and movement of fish into each regulatory area.
	| 
	|
	|
  	*/ 
  	int i,j,iyr;
  	
  	// average F to scale log_rbar
  	dvector fbar   = ("{0.05,0.15,0.07,0.15,0.07,0.03,0.05,0.05}");	

  	// Base initial abundance of halibut on the 2012 apportionment and move on.
  	// apportionment from 2012.
  	dvector apportionment(1,narea);
  	apportionment = ("{0.020,0.132,0.125,0.377,0.140,0.064,0.036,0.106}");

  	// std for log_rbar_devs
  	double sigma_R  = 0.0;

  	// std for fishery wpue and wbar
  	double sigma_it = 0.0;
  	double sigma_wt = 0.0;

  	// std for survey wpue
  	double sigma_yt = 0.0;

  	random_number_generator rng(seed);

  	Pj = moveProb;
  	// Pj = identity_matrix(1,narea);
  	growthModel();
  	initializeModel();

  	
  	dvector            ne(1,narea);
  	dvector     	   be(1,narea);
  	dvector            re(1,narea);
  	dmatrix M(1,narea,1,narea);
  	M = identity_matrix(1,narea);
  	M = value(Pj);
  	COUT(M)
  	
  	double s = exp(-value(m));
  	be = value(bo) * apportionment;
  	ne = be / value(wbar);   // Need to figure out equilribium soln for wbar
  	re = value((be-s*(alpha(syr)*ne+rho*be))/wk);
  	//ne = (value(bo) / value(wbar))/narea;
  	// re = value((ne - exp(-m)*ne))*M;
  	//re = ne*value(1.-exp(-m));
  	COUT(N*(1-s)*P);
  	COUT((N-s*N)*P);


  	COUT(wbar);
  	for( i = 1; i <= 200; i++ )
  	{
  		ne = (s*ne)*M + re;
  		be = s*value((alpha(syr)*ne+rho*be))*M + value(wk*re);
		cout<<s*(alpha(syr)*ne+rho*be)*M+value(wk*re)<<endl;
		//cout<<be<<endl;
  	}
  	re = ne - (s*ne)*M;
  	COUT(sum(re));
  	COUT(apportionment);
  	COUT(be/sum(be));
  	COUT(bo);
  	COUT(sum(re));
  	COUT(ro);

  	// COUT(bo*apportionment);
  	// exit(1);
  	// // Calculate average recruitment to each area (log_rbar)
  	// dvector   cbar(1,narea);
  	// dvector   wbar(1,narea);
  	
  	// dvector     se(1,narea);
  	
  	// cbar = colsum(obs_ct)/(nyr-syr+1);
  	// be   = elem_div(cbar , fbar);
  	// se   = exp( -value(m) - fbar );
  	// wbar = elem_div(se*value(alpha(syr)) + value(wk)*(1.0-se),(1.0-rho*se));
  	// re   = elem_prod(elem_div(be,wbar),1.0-se);
  	
  	// Add some random normal deviates to initial recruitment.
  	dvector   rtmp(1,narea);
 	rtmp.fill_randn(rng);
  	log_rbar = log(re)+sigma_R*rtmp;

  	// Random normal deviates to log_rbar_devs
  	// Ensure each column has a mean 0. Note that trans is not overloaded for
  	// init_bounded_matrix object. PAIN in the A**S
  	dmatrix r_devs(1,narea,syr,nyr);
  	r_devs.fill_randn(rng);
  	for( j = 1; j <= narea; j++ )
  	{
  		r_devs(j) -= mean(r_devs(j));
  		for( i = syr; i <= nyr; i++ )
  		{
  			log_rbar_devs(i,j) = sigma_R * r_devs(j,i);
  		}
  	}
  	
	
	// Initialize state variables
	initialStates();
	nt(syr) = ne;
	bt(syr) = be;

  	// Population dynamics
  	dvector v_bt(1,narea);
  	dvector v_ct(1,narea);
  	dvar_vector sj(1,narea);
	for(i=syr;i<=nyr;i++)
	{
		v_bt    = value(bt(i));
		v_ct    = 0*obs_ct(i);
		ft(i)   = getFt(value(m),v_bt,v_ct);
		sj      = mfexp( -m - ft(i) );
		nt(i+1) = elem_prod( sj, nt(i) )*Pj + rt(i);
		bt(i+1) = elem_prod( sj, value(alpha(i))*nt(i)+rho*bt(i) )*Pj + wk*rt(i);
	}
	COUT(rowsum(rt));
	// COUT(elem_div(bt,nt));
	COUT(bt);
	COUT(bt(nyr)/sum(bt(nyr)));
	exit(1);

	// Calculate fisheries catch statisitics and fill observations
	// with iid errors (it, and wt).
	calcFisheryCatchStatistics();
	dmatrix it_devs(1,n_fish_obs,1,narea);
	dmatrix wt_devs(1,n_wbar_obs,1,narea);
	it_devs.fill_randn(rng);
	wt_devs.fill_randn(rng);

	it_devs = sigma_it * it_devs;
	wt_devs = sigma_wt * wt_devs;

	// Overwrite WPUE data
	for( i = 1; i <= n_fish_obs; i++ )
	{
		iyr = fishery_wpue(i,0);
		for( j = 1; j <= narea; j++ )
		{
			double tmp = fishery_wpue(i,j);
			if( tmp > 0 )
			{
				fishery_wpue(i,j) = value(it(iyr,j))*exp(it_devs(i,j));
			}
		}
	}
	
	// Overwrite wbar data
	for( i = 1; i <= n_wbar_obs; i++ )
	{
		iyr = fishery_wbar(i,0);
		for( j = 1; j <= narea; j++ )
		{
			double tmp = fishery_wbar(i,j);
			if( tmp > 0 )
			{
				fishery_wbar(i,j) = value(wt(iyr,j))*exp(wt_devs(i,j));
			}
		}
	}

	// Calculate survey WPUE index
	calcSurveyStatistics();
	dmatrix yt_devs(1,n_surv_obs,1,narea);
	yt_devs.fill_randn(rng);
	yt_devs = sigma_yt * yt_devs;

	// Overwrite survey_wpue data
	for( i = 1; i <= n_surv_obs; i++ )
	{
		iyr = survey_wpue(i,0);
		for( j = 1; j <= narea; j++ )
		{
			double tmp = survey_wpue(i,j);
			if( tmp > 0 )
			{
				survey_wpue(i,j) = value(yt(iyr,j)) * exp(yt_devs(i,j));
			}
		}
	}
  }
	
FUNCTION calcFishingMortality
  {
	/*
		| This routine fills the matrix of fishing mortality rates
		| in each area and ensues each column of f_devs sums to zero.
	*/
	
	int i,j;
	
	for(i=syr;i<=nyr;i++)
	{
		ft(i) = mfexp( log_fbar + log_f_devs(i) );
	}
	
	// ft.sub(2000,2009) = 0.4;
	// ft(1996,4) = 0.5;
	
	/* penalty for log_f_devs */
	dvariable s;
	fpen.initialize();
	for(j=1;j<=narea;j++)
	{
		s     = mean(column(log_f_devs,j));
		fpen += 10000. * s*s;
	}
  }

FUNCTION growthModel
	/*
	| ----------------------------------------------------------------------- |
	| Calculate vector of alpha values based on a random walk                 |
	| ----------------------------------------------------------------------- |
	|
	*/
	int i;
	alpha = d_alpha;
	if( active(log_alpha) )
	{
		alpha = mfexp(log_alpha);
	}

	if( active(log_alpha_dev) )
	{
		for( i = syr; i < nyr; i++ )
		{
			alpha(i+1) = alpha(i)*exp(log_alpha_dev(i));
		}
	}
		
FUNCTION initializeModel
	dvariable s;
	
	// | Initialize coast-wide unfished states
	s    = mfexp(-m);
	bo   = mfexp(log_bo);
	reck = 4.0* steepness / (1.0 - steepness);
	wk   = alpha(syr)*( 1.0 - pow(rho,agek) )/(1.0-rho);
	wbar = (s*alpha(syr) + wk*(1.0-s))/(1.0-rho*s);
	no   = bo/wbar;
	ro   = no*(1.0-s);
	
	// | Parameters for stock-recruitment relationship
	// | Beverton holt model R=so*S/(1+beta*S)
	so   = reck*ro/bo;
	beta = (reck-1.0)/bo;
	

FUNCTION initialStates
	int i,j;
	// | Initialize recruitment vectors in all years and areas
	rbar = mfexp(log_rbar);
	for(i=syr; i<=nyr; i++)
	{
		rt(i) = mfexp(log_rbar + log_rbar_devs(i));
	}
	/* penalty for log_rbar_devs */
	dvariable s2;
	rpen.initialize();
	for(j=1;j<=narea;j++)
	{
		s2     = mean(column(log_rbar_devs,j));
		rpen += 1000. * s2*s2;
	}

	
	// | Initialize state variables for each regulatory area.
	dvar_vector sj(1,narea);
	dvar_vector wj(1,narea);
	sj      = mfexp(-m - ft(syr));  //FIXME should be fe
	wj      = elem_div( sj*alpha(syr) + wk*(1.0-sj),(1.0-rho*sj) );
	nt(syr) = elem_div( rt(syr),(1.0-sj) );
	bt(syr) = elem_prod(nt(syr),wj);

	// iteration for initialization to allow movement to stabilize.
	// for(int iter = 1; iter <= 50; iter++ )
	// {
	// 	nt(syr) = elem_prod(nt(syr),sj)*Pj + rt(syr);
	// 	bt(syr) = elem_prod( sj, alpha(syr)*nt(syr)+rho*bt(syr) )*Pj + wk*rt(syr);
	// }

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
	// p           = recruitmentDispersal;
	// COUT(recruitmentDispersal);
	
FUNCTION calcGravityModel
	/*
	|--------------------------------------------------------------------------|
	| MOVEMENT MODEL BASED ON GRAVITY WEIGTHS AND RESIDENCY 
	| This is referred to as a Markov transition matrix.               
	| -reference is Carruthers et al. 2011 and Caddy 1975                  
	| -estimate a vector of area-specific gravity weights and a residencey 
	|  parameter (probability of staying put in all areas).                
	| - set M_{f,f'} = r+g_{f}  if f = f'                                  
	|   set M_{f,f'} = g_{f}    if f!= f'                                  
	|       M_{f,f'} = M_{f,f'}/sum_{f'} M_{f,f'}                          
	| - M_{f,f'} is the probability of movement from f to f'
	|--------------------------------------------------------------------------|   
	*/
	int j;
	Pj.initialize();

	Delta   = (log_residency);
	gravity = (log_gravity_weights);

	// Delta = 5.0;
	// gravity.fill_seqadd(-2.0,0.4);
	dvar_matrix  G(1,narea,1,narea);
	dvar_matrix tG(1,narea,1,narea);

	for(j = 1; j <= narea; j++)
	{
		G(j,j) = mfexp(Delta + gravity(j));

		if( j < narea )
			G(j)(j+1,narea) = mfexp(gravity(j+1,narea));

		if( j > 1 )
			G(j)(1,j-1)     = mfexp(gravity(1,j-1));
	}
	
	tG = trans(G);
	dvar_vector colsumG=colsum(tG);
	for(j=1;j<=narea;j++)
	{
		G(j)=G(j)/colsumG(j);
	}
	Pj = (G);

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
		nt(i+1) = elem_prod( sj, nt(i) )*Pj + rt(i);
		bt(i+1) = elem_prod( sj, alpha(i)*nt(i)+rho*bt(i) )*Pj + wk*rt(i);
	}
  }
	
FUNCTION calcFisheryCatchStatistics
  {
	/*
	|-------------------------------------------------------------------------|
	| Compute the predicted fishery catch statistics:                         |
	|-------------------------------------------------------------------------|
	| 1) commercial catch by area          (ct)
	| 2) commercial wpue by area           (it)
	| 3) mean weight of the catch by area  (wt)
	*/
	
	int i,j,iyr;
	dvar_vector z(1,narea);
	dvar_vector u(1,narea);
	
	q = mfexp( log_q + log_q_dev - log(habitatArea) );
	for(i=syr;i<=nyr;i++)
	{
		/*1) commercial catch by area */
		z             = m + ft(i);
		u             = elem_prod(elem_div(ft(i),z),1.-exp(-z));
		ct(i)         = elem_prod(u,bt(i));
		varphi(i) = log(obs_ct(i)+TINY) - log(ct(i)+TINY);

		/*2) commercial wpue */
		// it(i) = elem_prod(q,bt(i));
		// or
		it(i)  = elem_div(elem_prod(q,ct(i)),ft(i));
		
		/*3) mean weight of the catch
			wt = (catch weight)/(catch numbers)
		*/
		// wt(i) = elem_div( ct(i),elem_prod(u,nt(i)) );
		wt(i) = elem_div(bt(i),nt(i));
	}

	// Residuals for commercial wpue.
	eta.initialize();
	for( i = 1; i <= n_fish_obs; i++)
	{
		iyr = fishery_wpue(i,0);
		for(j = 1; j<= narea; j++)
		{
			double tmp = fishery_wpue(i,j);
			if( tmp > 0 )
			{
				eta(i,j) = log(fishery_wpue(i,j)) - log(it(iyr,j));
			}
		}
	}

	// Residuals for average weight of commercial catch.
	nu.initialize();
	for( i = 1; i <= n_wbar_obs; i++)
	{
		iyr = fishery_wbar(i,0);
		for(j = 1; j<= narea; j++)
		{
			double tmp = fishery_wbar(i,j);
			if( tmp > 0 )
			{
				nu(i,j) = log(fishery_wbar(i,j)) - log(wt(iyr,j));
			}
		}
	}
  }

FUNCTION calcSurveyStatistics
  {
	/*
		| Compute predicted survey observations:
		| 1) relative abundance index for survey (yt)
	*/
	
	int i,j,iyr;
	
	qs = mfexp( log_qs + log_qs_dev - log(habitatArea) );
	for(i=syr;i<=nyr;i++)
	{
		/*1) relative abundance index for survey (yt) */
		yt(i)  = elem_prod(qs , bt(i));
	}

	// Residuals for setline survey
	epsilon.initialize();
	for( i = 1; i <= n_surv_obs; i++ )
	{
		iyr = survey_wpue(i,0);
		for( j = 1; j <= narea; j++ )
		{
			double tmp = survey_wpue(i,j);
			if( tmp > 0 )
			{
				epsilon(i,j) = log(survey_wpue(i,j)) - log(yt(iyr,j));
			}
		}
	}
  }

FUNCTION calcStockRecruit
  {
  	/*
  	|--------------------------------------------------------------------------|
	| Calculate stock recruitment relationship based on all areas combined.    |
	|--------------------------------------------------------------------------|
	| R_i = sum_j rt_{i,j}
	| S_i = sum_j bt_{i,j}
	|
  	*/
  	int i,j;
  	St = rowsum(bt.sub(syr,nyr));
  	Rt = rowsum(rt.sub(syr,nyr));

  	// Predicted recruits based on spawning biomass
  	for(i=syr; i<=nyr; i++)
  	{
  		Rhat(i+agek) = so * St(i) / (1.0 + beta*St(i) );
  	}
  	delta(syr+agek,nyr) = log(Rt(syr+agek,nyr)) - log(Rhat(syr+agek,nyr));
  }


FUNCTION calcObjectiveFunction
  {
  	int i,j;
  	/*
		Penalties to regularize the solution.
			1) pvec(1) = average fishing mortality ~ 0.2
			2) pvec(2) = norm2 on log_rec_devs

  	*/
  	dvar_vector	pvec(1,6);
  	pvec.initialize();
  	if( !last_phase() )
  	{
  		for( j=1; j<=narea; j++)
  		{
  			dvariable fbar = mean(column(ft,j));
  			pvec(1) += dnorm(log(fbar),log(0.2),0.10);
  			pvec(2) += 100. * norm2(column(log_f_devs,j));
  			pvec(3) += 100. * norm2(column(log_rbar_devs,j));
  		}
		pvec(4) = 100. * norm2(log_q_dev);
		pvec(5) = 100. * norm2(log_qs_dev);
		if( active(log_alpha_dev) )
		{
			pvec(6) = 100. * norm2(log_alpha_dev);
		}
  	}
  	else 
  	{
  		for( j=1; j<=narea; j++)
  		{
  			dvariable fbar = mean(column(ft,j));
  			pvec(1) += dnorm(log(fbar),log(0.2),1.0);
  			pvec(2) += 1.0 * norm2(column(log_f_devs,j));
  			pvec(3) += 1.0 * norm2(column(log_rbar_devs,j));
  		}	
  		pvec(4) = 10.0  * norm2(log_q_dev);
		pvec(5) = 10.0  * norm2(log_qs_dev);
		if( active(log_alpha_dev) )
		{
			pvec(6) = 100.  * norm2(log_alpha_dev);
		}
  	}

  	/*
	| Priors
	| 1) prior for steepness
  	*/
  	dvar_vector prior(1,1);
  	prior.initialize();

  	prior(1) = dbeta((steepness-0.2)/0.8,1.01,1.01);


  	/*
	| Negative loglikelihoods:
	| 	1) ell_catch = likelihood for the observed catch.
	| 	2) ell_wpue_fsh = likelihood for commercial cpue.
	| 	3) ell_wbar_fsh = likelihood for commercial wbar.
	| 	4) ell_wpue_srv = likelihood for the setline survey.
	| 	5) ell_rec      = likelihood for recruitment relationship.
	| 	6) ell_Pj       = likelihood for movement matrix.
  	*/
	ell_catch.initialize();
	ell_wpue_fsh.initialize();
	ell_wpue_srv.initialize();
	ell_Pj.initialize();
	

	for ( j = 1; j <= narea; j++)
	{
		// Catch likelihood
		dvar_vector res1 = column(varphi,j);
		ell_catch(j)     = dnorm(res1,0.05);

		// Commercial wpue likelihood
		dvar_vector res2 = column(eta,j);
		ell_wpue_fsh(j)  = dnorm(res2,0.40);

		// Commercial wbar likelihood
		dvar_vector res3 = column(nu,j);
		ell_wbar_fsh(j)  = dnorm(res3,0.30);

		// Setline survey likelihood (wpue)
		dvar_vector res4 = column(epsilon,j);
		ell_wpue_srv(j)  = dnorm(res4,0.20);

		// Movement probability
		ell_Pj(j) = 100.* norm2(moveProb(j)-Pj(j));
	}

	dvariable ell_rec;
	ell_rec = dnorm(delta,0.2);
	
	f  = fpen + rpen;
	f += sum(pvec);
	f += sum(prior);
	f += sum(ell_catch);
	f += sum(ell_wpue_fsh);
	f += sum(ell_wbar_fsh);
	f += sum(ell_wpue_srv);
	f += sum(ell_Pj);
	f += ell_rec;

  }

REPORT_SECTION
	ivector iyr(syr,nyr);
	ivector iyrs(syr,nyr+1);
	iyr.fill_seqadd(syr,1);
	iyrs.fill_seqadd(syr,1);
	REPORT(iyr);
	REPORT(iyrs);
	REPORT(agek);
	REPORT(alpha);
	REPORT(so);
	REPORT(beta);
	REPORT(St);
	REPORT(Rt);
	REPORT(Rhat);
	REPORT(Pj);
	REPORT(nt);
	REPORT(bt);
	REPORT(ft);
	REPORT(survey_wpue);
	REPORT(qs);
	REPORT(yt);
	REPORT(epsilon);
	REPORT(obs_ct);
	REPORT(ct);
	REPORT(varphi);
	REPORT(fishery_wbar);
	REPORT(wt);
	REPORT(nu);
	REPORT(fishery_wpue);
	REPORT(q);
	REPORT(it);
	REPORT(eta);
	REPORT(log_rbar_devs);
	REPORT(log_f_devs);



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
	#undef  REPORT
	#define REPORT(object) report << #object "\n" << object << endl;
	
	#undef  COUT
	#define COUT(object) cout<<setprecision(3)<<fixed<<#object "\n"<<object<<endl;
	
	#undef  TINY
	#define TINY 1.e-20

	#include <iostream>
	#include <iomanip>
	using namespace std;
	
	#include <admodel.h>
	#include <time.h>
	#include <statsLib.h>
	#include <Baranov.cpp>
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

