#include <stdio.h>
#include <string.h>
#include "math.h"
#include <stdarg.h>
#include <stdlib.h>
#include "../h/nr.h"
#include "../h/nrutil.h"
#include "../h/crust.h"
#include "../h/ns.h"
#include "../h/timer.h"

Crust::Crust() {

	// initialize default parameters

	// Neutron star parameters
	this->mass = 1.6;
	this->radius = 11.2;
	
	// Set up the grid
	this->N = 100;
	this->output = 1;
	this->hardwireQ=1;	
	this->B=0.0;
	this->accr=0;
	this->Tc = 3e7;
	this->yt = 1e12;
	
	// Envelope model
	this->gpe=0;	
	this->use_my_envelope=0;
	
	// Other parameters that are used during time evolution
	this->Qimp=1.0;
	this->rhot=1e9;
	this->rhob=1e14;	
	this->Tt=4e8;
	this->outburst_duration=(1.0/24.0) * 1.0/(365.0);
	this->extra_heating=0;
	this->force_cooling_bc=0;
	this->extra_Q=1.0;
	this->extra_y=1e12;
	this->deep_heating_factor=1.0;
	this->nuflag=1;
	this->Qrho=1e12;
	this->energy_slope=0.0;
	this->energy_deposited_inner=1.0;
	this->energy_deposited_outer=1.0;	
	this->timesofar=0.0;
	this->last_time_output=0.0;
	this->force_precalc=1;
	this->gap=1;
	this->kncrit=0.0;
	this->Lmin=0.0;
	this->Lscale=1.0;
	this->mdot=1.0;
	
	this->Qinner=-1.0;
	this->angle_mu=-1.0;
	
	this->use_potek_eos=0;
}

void Crust::setup(void) {

	if (this->yt < 10.0) this->yt=pow(10.0,this->yt);
	if (this->extra_y < 16.0) this->extra_y=pow(10.0,this->extra_y);

	if (this->Qinner == -1.0) this->Qinner=this->Qimp;
	if (this->energy_deposited_inner == -1.0) this->energy_deposited_inner = this->energy_deposited_outer;
	
	if (this->Qimp>=0.0) {   	// the Q values are assigned directly in 'calculate_vars'
		this->hardwireQ=1;
		printf("Using supplied Qimp values and HZ composition and heating.\n");
	} else {
		this->hardwireQ=0;
		printf("Using Qimp, composition, and heating from the crust model.\n");
	}

	if (this->angle_mu >= 0.0) this->B*=sqrt(0.75*this->angle_mu*this->angle_mu+0.25);
	printf("Magnetic field set to B=%lg\n", this->B);

	static Eos myeos(1);
	this->EOS = &myeos;
	this->EOS->Qimp=this->Qimp;
	this->EOS->gap=this->gap;
	this->EOS->kncrit=this->kncrit;
	this->EOS->B=this->B;
	this->EOS->accr=this->accr;
	this->EOS->use_potek_eos=this->use_potek_eos;
	
	set_ns_parameters(this->mass,this->radius,&this->g,&this->ZZ);
	set_up_grid("data/crust_model_shell");
	get_TbTeff_relation();
	
	// initialize the integrator
  	this->ODE.init(this->N+1,dynamic_cast<Ode_Int_Delegate *>(this));
	this->ODE.verbose=0;
  	this->ODE.stiff=1; this->ODE.tri=1;  // stiff integrator with tridiagonal solver
}

void Crust::evolve(double time, double mdot) {

	// time is the time to evolve in days
	// mdot is the accretion rate in Eddington units
	printf("Now evolve in time for %lg days at mdot=%lg\n",time,mdot);
	this->outburst_duration=time/(365.0*this->ZZ);
	if (mdot > 0.0) this->accreting = 1; else this->accreting = 0;
	this->mdot = mdot;

	clock_t timer;
	start_timing(&timer);
	int store_output=this->output;
	this->output=1;
  	precalculate_vars();
	this->output=store_output;
	this->force_precalc=0;  // only precalc once per session
	stop_timing(&timer,"precalculate_vars");

	start_timing(&timer);
	this->ODE.dxsav=1e4;
	for (int i=1; i<=this->N+1; i++) {
		this->ODE.set_bc(i,this->T[i]);
	}
	this->ODE.go(0.0, this->outburst_duration*3.15e7, this->outburst_duration*3.15e7*0.01,1e-6);
	stop_timing(&timer,"this->ODE.go");
	printf("number of steps = %d\n", this->ODE.kount);

	for (int i=1; i<=this->N+1; i++) {
		this->T[i]=ODE.get_y(i,this->ODE.kount);
	}

	// output results
	if (this->output) {
		printf("Starting output\n");
		if (this->last_time_output == 0.0) {
			this->fp=fopen("gon_out/out","w");
		   	this->fp2=fopen("gon_out/prof","w");
		} else {
			this->fp=fopen("gon_out/out","a");
	   		this->fp2=fopen("gon_out/prof","a");
		}
		if (this->last_time_output == 0.0) fprintf(this->fp,"%d %lg\n",this->N+1,this->g);
		start_timing(&timer);
		for (int j=1; j<=this->ODE.kount; j++) output_result_for_step(j,this->fp,this->fp2,this->timesofar,&this->last_time_output);
		fflush(this->fp); fflush(this->fp2);
		stop_timing(&timer,"output");
		fclose(this->fp);
		fclose(this->fp2);
		this->timesofar+=this->outburst_duration*3.15e7;
	}
}


Crust::~Crust() {
	// destructor
	this->ODE.tidy(); 
  	free_vector(this->rho,0,this->N+2);
  	free_vector(this->CP,0,this->N+2);
  	free_vector(this->P,0,this->N+2);
  	free_vector(this->K,0,this->N+2);
  	free_vector(this->F,0,this->N+2);
  	free_vector(this->NU,0,this->N+2);
  	free_vector(this->EPS,0,this->N+2);
}



void Crust::output_result_for_step(int j, FILE *fp, FILE *fp2,double timesofar,double *last_time_output) 
{
	// Output if enough time has elapsed
	if ((fabs(log10(fabs(timesofar+this->ODE.get_x(j))*this->ZZ)-log10(fabs(*last_time_output))) >= 0.01) || 
			(fabs(timesofar)+this->ODE.get_x(j))*this->ZZ < 1e5) {

		// get CP,K,eps,eps_nu at each point on the grid
		for (int i=1; i<=this->N+1; i++) calculate_vars(i,this->ODE.get_y(i,j),this->P[i],&this->CP[i],&this->K[i],&this->NU[i],&this->EPS[i]);

		// outer boundary
		double T0;
		outer_boundary(this->ODE.get_y(1,j),this->K[1],this->CP[1],this->NU[1],this->EPS[1],&T0,&this->K[0],&this->CP[0],&this->NU[0],&this->EPS[0]);

		// timestep
		double dt;
		if (j==1) dt=this->ODE.get_x(j); else dt=this->ODE.get_x(j)-this->ODE.get_x(j-1);

		// heat fluxes on the grid
		double *TT;
		TT=vector(1,this->N+1);
		for (int i=1; i<=this->N+1; i++) TT[i]=this->ODE.get_y(i,j);
		double FF = calculate_heat_flux(1,TT);
		for (int i=1; i<=this->N+1; i++) this->F[i] = calculate_heat_flux(i,TT);
		free_vector(TT,1,this->N+1);

		// total neutrino luminosity
		double Lnu=0.0;
		for (int i=1; i<=this->N; i++) Lnu += this->NU[i]*this->dx*this->P[i]/this->g;
	
		// we output time, fluxes and TEFF that are already redshifted into the observer frame
		// gon_out/prof
		fprintf(fp2, "%lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg\n", (timesofar+this->ODE.get_x(j))*this->ZZ, 
			pow((this->radius/11.2),2.0)*this->F[2]/(this->ZZ*this->ZZ), pow((this->radius/11.2),2.0)*FF/(this->ZZ*this->ZZ),
			this->ODE.get_y(this->N-5,j), pow((this->g/2.28e14)*TEFF.get(this->ODE.get_y(1,j))/5.67e-5,0.25)/this->ZZ, 
			this->ODE.get_y(1,j), pow((this->g/2.28e14)*TEFF.get(this->ODE.get_y(1,j))/5.67e-5,0.25),
			pow((this->radius/11.2),2.0)*this->F[this->N+1]/(this->ZZ*this->ZZ),pow((this->radius/11.2),2.0)*this->F[this->N]/(this->ZZ*this->ZZ),
			4.0*PI*pow(1e5*this->radius,2.0)*Lnu/(this->ZZ*this->ZZ), dt);
			
		if ((fabs(log10(fabs(timesofar+this->ODE.get_x(j))*this->ZZ)-log10(fabs(*last_time_output))) >= 1000.0) ||
			(fabs(timesofar)+this->ODE.get_x(j))*this->ZZ < 1e10) {
			// temperature profile into gon_out/out
			fprintf(fp,"%lg\n",this->ZZ*(timesofar+this->ODE.get_x(j)));
			for (int i=1; i<=this->N+1; i++)
				fprintf(fp, "%lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg\n", 
					this->P[i], this->ODE.get_y(i,j), this->F[i], this->NU[i], this->g*(this->F[i+1]-this->F[i])/(this->dx*this->P[i]), this->rho[i], this->CP[i]*this->rho[i], 
					this->ODE.get_d(i,j),1e8*pow(this->P[i]/2.521967e17,0.25), this->K[i], 2.521967e-15*pow(this->ODE.get_y(i,j),4)/this->P[i],
					this->NU[i],this->EPS[i]);
 		}
		*last_time_output=(timesofar+this->ODE.get_x(j))*this->ZZ;

	}


}




void Crust::set_up_grid(const char *fname)
// allocates storage for the grid and also computes the density
// and composition at each grid point
{
	// number of grid points this->N has already been set
	this->Pb=6.5e32; // pressure at the crust/core boundary
	this->Pt=this->yt*2.28e14;   // pressure at the top   // note I need to use 2.28 here to get the correct match to the envelope

	Spline QiSpline;
	Spline QhSpline;
	if (!this->hardwireQ) {   // only need this if we read in the crust model
		FILE *fp=fopen(fname,"r");
	
		int npoints;
		fscanf(fp,"%d",&npoints);
		npoints--;
		printf("Crust model has %d points\n",npoints);
		
		double *Qi,*Qh,*AA,*ZZ,*P,*Yn;
		Qi=vector(1,npoints);
		Qh=vector(1,npoints);
		AA=vector(1,npoints);
		ZZ=vector(1,npoints);
		Yn=vector(1,npoints);
		P=vector(1,npoints);

		for (int i=1; i<=npoints; i++) {
			double dummy;
			fscanf(fp, "%lg %lg %lg %lg %lg %lg %lg %lg\n",
				&P[i],&dummy,&dummy,&Qh[i],&ZZ[i],&AA[i],&Qi[i],&Yn[i]);
			//printf("%lg %lg %lg %lg %lg\n", rho[i], Qh[i],ZZ[i],AA[i],Qi[i]);
			P[i]=log10(P[i]);
		}

		YnSpline.minit(P,Yn,npoints);
		ZZSpline.minit(P,ZZ,npoints);
		AASpline.minit(P,AA,npoints);
		QiSpline.minit(P,Qi,npoints);
		QhSpline.minit(P,Qh,npoints);

		free_vector(Qi,1,npoints);
		free_vector(Qh,1,npoints);
		free_vector(AA,1,npoints);
		free_vector(ZZ,1,npoints);		
		free_vector(Yn,1,npoints);		
		free_vector(P,1,npoints);		
		fclose(fp);

	}

  	// storage
  	this->rho=vector(0,this->N+2);  
  	this->CP=vector(0,this->N+2);
  	this->P=vector(0,this->N+2);
  	this->K=vector(0,this->N+2);
  	this->F=vector(0,this->N+2);
  	this->T=vector(0,this->N+2);
  	this->NU=vector(0,this->N+2);
  	this->EPS=vector(0,this->N+2);
  	this->Qheat=vector(0,this->N+2);
  	this->Qimpur=vector(0,this->N+2);

  	this->dx=log(this->Pb/this->Pt)/(this->N-1);   // the grid is equally spaced in log column
  
	FILE *fp=NULL;
	if (this->output) fp = fopen("gon_out/grid_profile","w");

	double Qtot=0.0;
  	for (int i=0; i<=this->N+2; i++) {
    	double x=log(this->Pt)+this->dx*(i-1);
    	this->P[i]=exp(x);
		this->EOS->P = this->P[i];
		  // we have to set the temperature to something
		this->T[i] = this->Tc;
		this->EOS->T8=this->T[i]/1e8; 
		set_composition();
		this->EOS->rho=this->EOS->find_rho();
		this->rho[i]=this->EOS->rho;

		// GammaT[i] refers to i+1/2
		// The following line uses a composition of 56Fe to calculate gamma,
		// it avoids jumps in the melting point associated with e-capture boundaries
		double GammaT;
		if (0) {
			double Z=26.0,A=56.0;   // 28Si
			GammaT = pow(Z*4.8023e-10,2.0)*pow(4.0*PI*this->EOS->rho/(3.0*A*1.67e-24),1.0/3.0)/1.38e-16;
		} else {
			GammaT = pow(this->EOS->Z[1]*4.8023e-10,2.0)*pow(4.0*PI*this->EOS->rho/(3.0*this->EOS->A[1]*1.67e-24),1.0/3.0)/1.38e-16;
		}

		double Tmelt = 5e8*pow(this->P[i]/(2.28e14*1.9e13),0.25)*pow(this->EOS->Z[1]/30.0,5.0/3.0);
		double LoverT = 0.8 * 1.38e-16 /(this->EOS->A[1]*1.67e-24);

		this->Qheat[i]=0.0;
		if (!this->hardwireQ) {
			double P1 = exp(x-0.5*this->dx);
			double P2 = exp(x+0.5*this->dx);
			this->Qheat[i]=QhSpline.get(log10(P2))-QhSpline.get(log10(P1));
			if (this->Qheat[i]<0.0) this->Qheat[i]=0.0;
		}
		Qtot+=this->Qheat[i];

		if (!this->hardwireQ) {
			this->Qimpur[i]=QiSpline.get(log10(this->P[i]));
			if (this->Qimpur[i] < 0.0) this->Qimpur[i]=0.0;
		}

		if (this->output) 
			fprintf(fp, "%d %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg\n", i, this->P[i], this->rho[i], this->EOS->A[1]*(1.0-this->EOS->Yn), 
				this->EOS->Z[1], this->EOS->Yn,this->EOS->A[1],this->EOS->ptot(), Tmelt, GammaT/1e8, LoverT*1e8);
			//,AASpline.get(log10(this->rho[i])),  ZZSpline.get(log10(this->rho[i])), this->Qimpur[i], this->Qheat[i]);
  	}

	if (this->output) fclose(fp);

	if (!this->hardwireQ) QiSpline.tidy();

  	printf("Grid has %d points, delx=%lg, Pb=%lg, rhob=%lg, Pt=%lg, rhot=%lg\n", 
			this->N, this->dx, this->P[this->N],this->rho[this->N],this->P[1],this->rho[1]);
		printf("Total heat release is %lg MeV\n",Qtot);
}



void Crust::get_TbTeff_relation(void)
// reads in the Flux-T relation from the data file output by makegrid.cc
{
	double *temp, *flux;  // temporary storage to initialize the spline
	int npoints = 195;  //  needs to be >= number of points read in
	temp = vector(1,npoints);
	flux = vector(1,npoints);
	
	// the file "out/grid" is made by makegrid.cc
	// it contains  (column depth, T, flux)  in cgs
	FILE *fp;
	if (this->use_my_envelope) {
		if (this->EOS->B == 1e15) fp = fopen("out/grid_1e15_nopotek","r");
		else if (this->EOS->B == 1e14) fp = fopen("out/grid_1e14_potek","r");
		else if (this->EOS->B == 3e14) fp = fopen("out/grid_3e14_potek","r");
		else if (this->EOS->B == 3e15) fp = fopen("out/grid_3e15_potek","r");
		else {
			printf("Don't know which envelope model to use for this B!\n");
			exit(1);
		}
	} else {
		if (this->gpe) fp = fopen("out/grid_He4","r");
		else fp = fopen("out/grid_He9","r");
	}
	FILE *fp2=NULL;
	if (this->output) fp2=fopen("gon_out/TbTeff", "w");
	
	double y,T,F,rho,dummy;
	int count = 0;
	while (!feof(fp)) {
		fscanf(fp, "%lg %lg %lg %lg %lg %lg\n", &y, &T, &F,&rho,&dummy,&dummy);
		if (fabs(y-log10(this->yt))<1e-3) {  // select out the points which correspond to the top column
			count++;
			temp[count] = pow(10.0,T);
			// correct for gravity here:
			flux[count] = pow(10.0,F); //* (this->g/2.28e14);
			if (this->output) fprintf(fp2, "%d %lg %lg %lg %lg %lg\n", count,y,T,F,temp[count],flux[count]);
		}
	}
		
	fclose(fp);
	if (this->output) fclose(fp2);
	
	// the following spline contains the flux as a function of temperature at column depth this->yt
	this->TEFF.minit(temp,flux,count);
	
	free_vector(temp,1,npoints);
	free_vector(flux,1,npoints);
}



void Crust::set_composition(void)
{
	if (this->hardwireQ) {
		// use the EOS routines to get the composition
		// ie. crust models from the literature
		this->EOS->set_composition_by_pressure();	
	} else {
		// otherwise use our crust model
		// the model gives the mean A, mean Z and Yn
		// and we set up the variables as in our this->EOS->set_comp() routine
		this->EOS->Yn = YnSpline.get(log10(this->EOS->P));
		if (this->EOS->Yn < 1e-6) this->EOS->Yn=0.0;
		this->EOS->A[1]=AASpline.get(log10(this->EOS->P));
		this->EOS->A[1]/=(1.0-this->EOS->Yn);
		this->EOS->Z[1]=ZZSpline.get(log10(this->EOS->P));
		this->EOS->set_Ye=this->EOS->Z[1]/this->EOS->A[1];		
		//printf("%lg %lg %lg %lg\n", this->EOS->Yn, this->EOS->rho, this->EOS->A[1], this->EOS->Z[1]);
	}
}



void Crust::precalculate_vars(void) 
// calculate various quantities at each grid point as a function of temperature
// then during the run we can look them up in a table
{
	// the table is constructed in terms of log10(T)
	// for historical reasons, this is called beta here
	// (for long X-ray bursts where radiation pressure is significant,
	// beta=Prad/P is a better variable to use)
	this->nbeta=100;
	this->betamin=6.5;
	this->betamax=10.0;
	this->deltabeta = (this->betamax-this->betamin)/(1.0*(this->nbeta-1));	

	this->CP_grid = matrix(1,this->N+2,1,this->nbeta);	
	this->K1_grid = matrix(1,this->N+2,1,this->nbeta);	
	this->K0_grid = matrix(1,this->N+2,1,this->nbeta);	
	this->KAPPA_grid = matrix(1,this->N+2,1,this->nbeta);	
	this->K1perp_grid = matrix(1,this->N+2,1,this->nbeta);	
	this->K0perp_grid = matrix(1,this->N+2,1,this->nbeta);	
	this->NU_grid = matrix(1,this->N+2,1,this->nbeta);	
	this->EPS_grid = matrix(1,this->N+2,1,this->nbeta);	

	// For the crust heating, we need to convert the density limits into 
	// pressures
	EOS->rho = this->rhot;
	EOS->set_composition_by_density();
	this->heating_P1 = EOS->ptot();
	EOS->rho = this->rhob;
	EOS->set_composition_by_density();
	this->heating_P2 = EOS->ptot();

	FILE *fp = NULL;
	char s[100];
	if (EOS->B > 0.0) sprintf(s,"gon_out/precalc_results_%lg",log10(EOS->B));
	else sprintf(s,"gon_out/precalc_results_0");
	if (!this->force_precalc) fp=fopen(s,"r");
	// if unsuccessful (or if precalc is set) we need to recalculate
	if (fp == NULL) {
		if (this->output)
			fp=fopen(s,"w");
		printf("Precalculating quantities and writing to file %s...\n",s);

		for (int i=1; i<=this->N+1; i++) {
	
			EOS->P=this->P[i];
			EOS->rho = this->rho[i];
			set_composition();
			
			if (this->output) 
				fprintf(fp, "Grid point %d  P=%lg  rho=%lg  A=%lg  Z=%lg Yn=%lg:  T8,CP,K,eps_nu,eps_nuc\n",
					i, this->P[i], this->rho[i], (1.0-EOS->Yn)*EOS->A[1], EOS->Z[1], EOS->Yn);
		
			double heating = crust_heating(i);
		
			for (int j=1; j<=this->nbeta; j++) {		
				double beta = this->betamin + (j-1)*(this->betamax-this->betamin)/(1.0*(this->nbeta-1));
				EOS->T8 = 1e-8*pow(10.0,beta);

				this->CP_grid[i][j]=EOS->CV();
				this->NU_grid[i][j]=EOS->eps_nu();			
				this->EPS_grid[i][j]=heating;

				// we calculate the thermal conductivity for Q=0 and Q=1, and later interpolate to the
				// current value of Q. This means we can keep the performance of table lookup even when
				// doing MCMC trials which vary Q.

				double Q_store=EOS->Qimp;  // store Q temporarily

				EOS->Qimp=0.0;
				double Kcond,Kcondperp;
				//Kcond = EOS->K_cond(EOS->Chabrier_EF());
				//Kcondperp=Kcond;
				Kcond = EOS->potek_cond();
				Kcondperp = EOS->Kperp;   
				this->K0_grid[i][j]=EOS->rho*Kcond/this->P[i];
				this->K0perp_grid[i][j]=EOS->rho*Kcondperp/this->P[i];

				EOS->Qimp=1.0;
				//Kcond = EOS->K_cond(EOS->Chabrier_EF());
				//Kcondperp=Kcond;
				Kcond = EOS->potek_cond();
				Kcondperp = EOS->Kperp;
				this->K1_grid[i][j]=EOS->rho*Kcond/this->P[i];
				this->K1perp_grid[i][j]=EOS->rho*Kcondperp/this->P[i];

				EOS->Qimp=Q_store;  // restore to previous value

				// conductivity due to radiation
				(void) EOS->opac();  // call to kappa sets the variable kappa_rad
				this->KAPPA_grid[i][j] = 3.03e20*pow(EOS->T8,3)/(EOS->kappa_rad*this->P[i]);

				if (this->output) 
					fprintf(fp, "%lg %lg %lg %lg %lg %lg %lg %lg %lg\n", EOS->T8, this->CP_grid[i][j], 
						this->K0_grid[i][j],this->K1_grid[i][j], this->K0perp_grid[i][j],this->K1perp_grid[i][j],
						this->NU_grid[i][j], this->EPS_grid[i][j], this->KAPPA_grid[i][j] );

				this->EPS_grid[i][j]*=energy_deposited(i);

			}	
		}
		if (this->output) fclose(fp);

	} else {
		
		printf("***Reading precalculated quantities from file %s\n", s);
		
		for (int i=1; i<=this->N+1; i++) {
			int kk; double dd;
			fscanf(fp, "Grid point %d  P=%lg  rho=%lg  A=%lg  Z=%lg Yn=%lg:  T8,CP,K,eps_nu,eps_nuc\n",
					&kk,&dd,&dd,&dd,&dd,&dd);
			for (int j=1; j<=this->nbeta; j++) {		
				fscanf(fp, "%lg %lg %lg %lg %lg %lg %lg %lg %lg\n", &EOS->T8, &this->CP_grid[i][j], 
					&this->K0_grid[i][j],&this->K1_grid[i][j], &this->K0perp_grid[i][j],&this->K1perp_grid[i][j],
					&this->NU_grid[i][j], &this->EPS_grid[i][j],&this->KAPPA_grid[i][j]);
				// always calculate the crust heating..
				this->EPS_grid[i][j]=crust_heating(i);
				this->EPS_grid[i][j]*=energy_deposited(i);
			}
		}
		fclose(fp);
			
	}
	
}


double Crust::energy_deposited(int i)
{
	double ener;
	if (this->rho[i]>4e11) ener = this->energy_deposited_inner;
	else ener = this->energy_deposited_outer;
	ener *= pow(this->rho[i]/1e10,this->energy_slope);
	return ener;
}


double Crust::crust_heating(int i) 
// calculates the crust heating in erg/g/s
// for grid point i
{
	double eps=0.0,P = this->P[i];

	// if we are heating on < 1day timescale then its a magnetar
	if (this->outburst_duration<1.0/365.0) {
		
		// eps in erg/g/s
		double eps_heat = 1e25/(this->rho[i]*this->outburst_duration*3.15e7);
		eps_heat /= this->mdot * this->g;   // modify to the units used in the code

		// limit the heating to a region of the crust
		double P1 = P*exp(-0.5*this->dx);
		double P2 = P*exp(0.5*this->dx);
		if (P1 > this->heating_P1 && P2 < this->heating_P2)   // we are within the heating zone
			eps = eps_heat;
		if (P1 < this->heating_P1 && P2 < this->heating_P2 && this->heating_P1<P2) {   // left hand edge of heated region
			eps = eps_heat * log(P2/this->heating_P1)/this->dx;
		}
		if (P1 > this->heating_P1 && P2 > this->heating_P2 && this->heating_P2>P1) {  // right hand edge of heated region
			eps = eps_heat * log(this->heating_P2/P1)/this->dx;
		}

	} else {  // otherwise we are doing an accreting neutron star

		if (!this->hardwireQ) {
			eps = this->Qheat[i]*8.8e4*9.64e17/(this->P[i]*this->dx);
		} else {
			// simple "smeared out" heating function, 1.2MeV in inner crust, 0.2MeV in outer crust
			if (P >= 1e16*2.28e14 && P <= 1e17*2.28e14) eps=8.8e4*this->deep_heating_factor*1.7*9.64e17/(P*log(1e17/1e16));
		 	if (P >= 3e12*2.28e14 && P < 3e15*2.28e14) eps=8.8e4*this->deep_heating_factor*0.2*9.64e17/(P*log(3e15/3e12));

			// Extra heat source in the ocean
			if (this->extra_heating) {	
				// Put all of the extra heat into one grid point
				//if (this->P[i]*exp(-0.5*this->dx) <this->extra_y*2.28e14 && this->P[i]*exp(0.5*this->dx)>this->extra_y*2.28e14)
				//		eps+=8.8e4*this->extra_Q*9.64e17/(P*this->dx);

				// More distributed heating
				double extra_y1 = this->extra_y/3.0;	
				double extra_y2 = this->extra_y*3.0;
				double eps_extra=0.0;
				double P1 = P*exp(-0.5*this->dx);
				double P2 = P*exp(0.5*this->dx);
				double geff=2.28e14;

				if (P1 > extra_y1*geff && P2 < extra_y2*geff)   // we are within the heating zone
					eps_extra=8.8e4*this->extra_Q*9.64e17/(P*log(extra_y2/extra_y1));
				if (P1 < extra_y1*geff && P2 < extra_y2*geff && extra_y1*geff<P2) {   // left hand edge of heated region
					eps_extra=8.8e4*this->extra_Q*9.64e17/(P*log(extra_y2/extra_y1));
					eps_extra *= log(P2/(extra_y1*geff))/this->dx;
				}
				if (P1 > extra_y1*geff && P2 > extra_y2*geff && extra_y2*geff>P1) {  // right hand edge of heated region
					eps_extra=8.8e4*this->extra_Q*9.64e17/(P*log(extra_y2/extra_y1));
					eps_extra *= log(extra_y2*geff/P1)/this->dx;
				}
				eps+=eps_extra;
			}
		}
	}

	return eps;	
}





// ----------------------------- Time integration -------------------------------------------


void Crust::derivs(double t, double T[], double dTdt[])
// calculates the time derivatives for the whole grid
{
	// First calculate quantities at each grid point
	for (int j=1; j<=this->N; j++) calculate_vars(j,T[j], this->P[j], &this->CP[j], &this->K[j], &this->NU[j],&this->EPS[j]);
  	outer_boundary(T[1],this->K[1],this->CP[1],this->NU[1],this->EPS[1],&T[0],&this->K[0],&this->CP[0],&this->NU[0],&this->EPS[0]);
  	inner_boundary(T[this->N],this->K[this->N],this->CP[this->N],this->NU[this->N],this->EPS[this->N],
				&T[this->N+1],&this->K[this->N+1],&this->CP[this->N+1],&this->NU[this->N+1],&this->EPS[this->N+1]);

	// determine the fluxes at the half-grid points
	//  this->F[i] is the flux at i-1/2
  	for (int i=1; i<=this->N+1; i++)   this->F[i] = calculate_heat_flux(i,T);	
	
	// Calculate the derivatives dT/dt
	for (int i=1; i<=this->N; i++) {
  		dTdt[i]=this->g*(this->F[i+1]-this->F[i])/(this->dx*this->CP[i]*this->P[i]);
		if (this->nuflag) dTdt[i]+=-(this->NU[i]/this->CP[i]);
		if (this->accreting) dTdt[i]+=this->EPS[i]/this->CP[i];
	}
  	dTdt[this->N+1]=0.0;
}

double Crust::calculate_heat_flux(int i, double *T)
{
	double flux;
	if (i>1 || (this->accreting && this->outburst_duration > 1.0/365.0 && !this->force_cooling_bc))
//		if (i>1 || (this->accreting && EOS->B == 0.0))   
		// use this inside the grid, or at the surface when we are accreting (which 
		// fixes the outer temperature)
	 	flux = 0.5*(this->K[i]+this->K[i-1])*(T[i]-T[i-1])/this->dx;	
	else {
		// cooling boundary condition
		if (EOS->B == 0.0 || this->use_my_envelope) {
			// from my envelope calculation (makegrid.cc)
			flux = (this->g/2.28e14)*TEFF.get(T[i]);
		} else {
			// for magnetars we use
			// Potekhin & Yakovlev 2001 eq.(27)
			double T9 = T[i]*1e-9;
			double xi = T9 - 0.001*pow(1e-14*this->g,0.25)*sqrt(7.0*T9);
			flux = 5.67e-5 * 1e24 * this->g*1e-14 * (pow(7*xi,2.25)+pow(0.333*xi,1.25));
		
			// or use makegrid.cc calculation
			//flux = (this->g/2.28e14)*TEFF.get(T[i]);
			
			
			// now correct for B ... 
			if (this->angle_mu >= 0.0) {
				// use the enhancement along the field direction
				double B12=EOS->B*1e-12;
				double chi1 = 1.0 + 0.0492*pow(B12,0.292)/pow(T9,0.24);
				//double chi2 = sqrt(1.0 + 0.1076*B12*pow(0.03+T9,-0.559))/
				//			pow(1.0+0.819*B12/(0.03+T9),0.6463);
				double fcond = 4.0*this->angle_mu*this->angle_mu/(1.0+3.0*this->angle_mu*this->angle_mu);		
				flux *= fcond*pow(chi1,4.0);//+(1.0-fcond)*pow(chi2,4.0);

			} else {
				// or use eq. (31) or PY2001  which gives F(B)/F(0)
				double fac, a1,a2,a3,beta;
				beta = 0.074*sqrt(1e-12*EOS->B)*pow(T9,-0.45);
				a1=5059.0*pow(T9,0.75)/sqrt(1.0 + 20.4*sqrt(T9) + 138.0*pow(T9,1.5) + 1102.0*T9*T9);
				a2=1484.0*pow(T9,0.75)/sqrt(1.0 + 90.0*pow(T9,1.5)+ 125.0*T9*T9);
				a3=5530.0*pow(T9,0.75)/sqrt(1.0 + 8.16*sqrt(T9) + 107.8*pow(T9,1.5)+ 560.0*T9*T9);
				fac = (1.0 + a1*beta*beta + a2*pow(beta,3.0) + 0.007*a3*pow(beta,4.0))/(1.0+a3*beta*beta);
				flux *= fac;
			}
		}
	}
		
	return flux;
}

double Crust::dTdt(int i, double *T)
// calculates the time derivative for grid cell i 
// This is used when calculating the jacobian in tri-diagonal form
{
	int k=i-1; if (k<1) k=1;
	int k2=i+1; if (k2>this->N+1) k2=this->N+1;
	for (int j=k; j<=k2; j++) calculate_vars(j,T[j], this->P[j], &this->CP[j], &this->K[j], &this->NU[j],&this->EPS[j]);
	if (i==1) outer_boundary(T[1],this->K[1],this->CP[1],this->NU[1],this->EPS[1],&T[0],&this->K[0],&this->CP[0],&this->NU[0], &this->EPS[0]);
	if (i==this->N) inner_boundary(T[this->N],this->K[this->N],this->CP[this->N],this->NU[this->N], this->EPS[this->N],
										&T[this->N+1],&this->K[this->N+1],&this->CP[this->N+1],&this->NU[this->N+1],&this->EPS[this->N+1]);

	double f=this->g*(calculate_heat_flux(i+1,T)-calculate_heat_flux(i,T))/(this->dx*this->CP[i]*this->P[i]);
	if (this->nuflag) f+=-(this->NU[i]/this->CP[i]);
	if (this->accreting) f+=this->EPS[i]/this->CP[i];

	return f;
}

void Crust::outer_boundary(double T1, double K1, double CP1, double NU1, double EPS1,
		double *T0, double *K0, double *CP0, double *NU0, double *EPS0)  
{
	if (this->accreting && this->Tt>0.0) *T0=this->Tt;   // constant temperature during accretion
	else *T0=T1*(8.0-this->dx)/(8.0+this->dx);   // assumes radiative zero solution, F\propto T^4
	*K0=K1; *CP0=CP1;
	if (this->nuflag) *NU0=NU1; else *NU0=0.0;
	if (this->accreting) *EPS0=EPS1; else *EPS0=0.0;
}

void Crust::inner_boundary(double TN, double KN, double CPN, double NUN, double EPSN,
	double *TN1, double *KN1, double *CPN1, double *NUN1, double *EPSN1)
{
	*TN1=this->Tc;   // fixed core temperature
//	*TN1=TN;    // zero flux inner boundary
	*KN1=KN;
	*CPN1=CPN;  
	if (this->nuflag) *NUN1=NUN; else *NUN1=0.0;
	if (this->accreting) *EPSN1=EPSN; else *NUN1=0.0;
	
}

void Crust::jacobn(double t, double *T, double *dfdt, double **dfdT, int n)
// calculates the Jacobian numerically
{
  	double e=0.01;

	// takes advantage of the tri-diagonal nature to calculate as few dTdt's as needed
	double f;
  // this assumes the arrays dfdt and dfdT are preinitialized to zero (I changed odeint to do this)
  for (int i=2; i<n; i++) {
  	T[i-1]*=1.0+e; f=dTdt(i,T);
    T[i-1]/=1.0+e; dfdT[i][i-1]=(f-dfdt[i])/(T[i-1]*e);
    T[i]*=1.0+e; f=dTdt(i,T);
    T[i]/=1.0+e; dfdT[i][i]=(f-dfdt[i])/(T[i]*e);
    T[i+1]*=1.0+e; f=dTdt(i,T);
    T[i+1]/=1.0+e; dfdT[i][i+1]=(f-dfdt[i])/(T[i+1]*e);
  }
  {
	int i=1;
	T[i]*=1.0+e; f=dTdt(i,T);
	T[i]/=1.0+e; dfdT[i][i]=(f-dfdt[i])/(T[i]*e);
	T[i+1]*=1.0+e; f=dTdt(i,T);
	T[i+1]/=1.0+e; dfdT[i][i+1]=(f-dfdt[i])/(T[i+1]*e);
  }
}  


void Crust::calculate_vars(int i, double T, double P, double *CP, double *K, double *NU,double *EPS)
{
	// sometimes we get a nan value for T here from the integrator
	// In this case, set the temperature to be some value.. this seems to
	// deal with this problem ok
	if (isnan(T) || T<0.0) T=1e7;
	
	double beta=log10(T);
	// if beta lies outside the table, set it to the max or min value
	if (beta > this->betamax) beta = this->betamax;
	if (beta < this->betamin) beta = this->betamin;
		
		// lookup values in the precalculated table
	int j = 1 + (int) ((beta-this->betamin)/this->deltabeta);
	double interpfac=(beta-(this->betamin + (j-1)*this->deltabeta))/this->deltabeta;
	// interpolate the thermal conductivity to the current
	// value of impurity parameter Q
	double K0=this->K0_grid[i][j] + (this->K0_grid[i][j+1]-this->K0_grid[i][j])*interpfac;
	double K1=this->K1_grid[i][j] + (this->K1_grid[i][j+1]-this->K1_grid[i][j])*interpfac;
	//double K0perp=this->K0perp_grid[i][j] + (this->K0perp_grid[i][j+1]-this->K0perp_grid[i][j])*interpfac;
	//double K1perp=this->K1perp_grid[i][j] + (this->K1perp_grid[i][j+1]-this->K1perp_grid[i][j])*interpfac;
	//K0perp=0.0; K1perp=0.0;

	// use something like this next line to hardwire Q values
	double Qval;
	if (this->hardwireQ) {
		if (this->rho[i] > this->Qrho) Qval=this->Qinner; else Qval=EOS->Qimp;
//		if (P>2.28e29) Qval=this->Qinner; else Qval=EOS->Qimp;
	} else {
		Qval = this->Qimpur[i];	
	}
	double KK,KKperp;
	KK=this->g*K0*K1/(K0*Qval+(1.0-Qval)*K1);

	double kappa;
	kappa=this->KAPPA_grid[i][j] + (this->KAPPA_grid[i][j+1]-this->KAPPA_grid[i][j])*interpfac;
	kappa*=this->g;
	KK += kappa;
	
	if (EOS->B > 0) {
		KKperp=0.0; //this->g*K0perp*K1perp/(K0perp*Qval+(1.0-Qval)*K1perp);
		if (this->angle_mu >= 0.0) {
			KK *= 4.0*this->angle_mu*this->angle_mu/(1.0+3.0*this->angle_mu*this->angle_mu);
		} else {
			KK = 0.5*(1.0544*KK+0.9456*KKperp);  // average over dipole geometry	
		}
	}
//	if (EOS->B > 0) {
//	//KKperp = this->g*K0perp*K1perp/(K0perp*Qval+(1.0-Qval)*K1perp);		
//	double fcond = 4.0*this->angle_mu*this->angle_mu/(1.0+3.0*this->angle_mu*this->angle_mu);		
//	*K=fcond*KK;//+(1.0-fcond)*KKperp;	
//} else {
	*K=KK;
//}
	
	*CP=this->CP_grid[i][j] + (this->CP_grid[i][j+1]-this->CP_grid[i][j])*interpfac;
	if (this->nuflag) *NU=this->NU_grid[i][j] + (this->NU_grid[i][j+1]-this->NU_grid[i][j])*interpfac; 
	else *NU=0.0;
	if (this->accreting) {
		*EPS=this->EPS_grid[i][1];  // assume heating is independent of temperature 
	//	*EPS=(this->EPS_grid[i][j] + (this->EPS_grid[i][j+1]-this->EPS_grid[i][j])*interpfac); 
		*EPS=*EPS * this->mdot * this->g;
	}
	else *EPS=0.0;
 }





