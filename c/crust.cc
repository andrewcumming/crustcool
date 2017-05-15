#include <stdio.h>
#include <string.h>
#include "math.h"
#include <stdarg.h>
#include <stdlib.h>
#include "../h/ns.h"
#include "../h/vector.h"
#include "../h/crust.h"
#include "../h/timer.h"

// --------------------------------- Constructor and destructor ---------------------------------------------

Crust::Crust() {

	// initialize the default parameters

	// Neutron star parameters
	this->mass = 1.6;
	this->radius = 11.2;
	
	// Heat capacity and neutrino emissivity of the core (innermost grid point)
	this->C_core = 1e38;
	this->Lnu_core_norm = 1e31;
	this->Lnu_core_alpha = 8.0;
		
	// Set up the grid
	this->N = 100;
	this->output = 1;
	this->hardwireQ=1;	// default =1 means that Qimp is specified in the init.dat file; otherwise Qimp(rho) is read in from the crust model
	this->B=0.0;
	this->accr=0;    // accreted or non-accreted crust?
	this->Tc = 3e7;
	this->yt = 1e12;
	
	// Envelope model
	this->gpe=0;	
	this->use_my_envelope=0;
	
	// Other parameters that are used during time evolution
	this->Qimp=1.0;
	this->rhot=1e9;
	this->Tt=4e8;
	this->rhob=1e14;
	this->outburst_duration=(1.0/24.0) * 1.0/(365.0);   // default outburst is 1 hour => magnetar
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
	
	this->resume = 0;    // if =1 then read in the temperature profile from last time and start from there
}


Crust::~Crust() {
	// destructor
	this->ODE.tidy(); 
	delete [] this->grid;
}

// --------------------------------- Setup ---------------------------------------------

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
		Qi=new double [npoints+1];
		Qh=new double [npoints+1];
		AA=new double [npoints+1];
		ZZ=new double [npoints+1];
		Yn=new double [npoints+1];
		P=new double [npoints+1];

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

		delete [] Qi;
		delete [] Qh;
		delete [] AA;
		delete [] ZZ;
		delete [] Yn;
		delete [] P;
		fclose(fp);

	}

  	// storage
	this->grid = new GridPoint [this->N+2];

	// grid spacing (equal spacing in log column)
 	this->dx=log(this->Pb/this->Pt)/(this->N-1);
  
	FILE *fp=NULL;
	if (this->output) fp = fopen("out/grid_profile","w");

	double Qtot=0.0;
  	for (int i=0; i<=this->N+1; i++) {
    	double x=log(this->Pt)+this->dx*(i-1);
    	this->grid[i].P=exp(x);
		this->EOS->P = this->grid[i].P;
		  // we have to set the temperature to something
		this->grid[i].T = this->Tc;
		this->EOS->T8=this->grid[i].T/1e8; 
		set_composition();
		this->EOS->rho=this->EOS->find_rho();
		this->grid[i].rho=this->EOS->rho;

		// integrate the equation of hydrostatic balance to get
		// the radial location of each grid point
		if (i==0) {
			this->grid[i].r = this->radius * 1e5;
		} else {
			this->grid[i].r = this->grid[i-1].r - this->dx * this->grid[i].P/(this->ZZ*this->grid[i].rho*this->g);
		}
	
		// GammaT[i] refers to i+1/2
		// The following line uses a composition of 56Fe to calculate gamma,
		// it avoids jumps in the melting point associated with e-capture boundaries
		double GammaT;
		if (0) {
			double Z=26.0,A=56.0;   // 56Fe
			GammaT = pow(Z*4.8023e-10,2.0)*pow(4.0*M_PI*this->EOS->rho/(3.0*A*1.67e-24),1.0/3.0)/1.38e-16;
		} else {
			GammaT = pow(this->EOS->Z[1]*4.8023e-10,2.0)*pow(4.0*M_PI*this->EOS->rho/(3.0*this->EOS->A[1]*1.67e-24),1.0/3.0)/1.38e-16;
		}

		double Tmelt = 5e8*pow(this->grid[i].P/(2.28e14*1.9e13),0.25)*pow(this->EOS->Z[1]/30.0,5.0/3.0);
		double LoverT = 0.8 * 1.38e-16 /(this->EOS->A[1]*1.67e-24);

		this->grid[i].Qheat=0.0;
		if (!this->hardwireQ) {    // we're using our own crust model
			double P1 = exp(x-0.5*this->dx);
			double P2 = exp(x+0.5*this->dx);
			this->grid[i].Qheat=QhSpline.get(log10(P2))-QhSpline.get(log10(P1));
			if (this->grid[i].Qheat<0.0) this->grid[i].Qheat=0.0;
			Qtot+=this->grid[i].Qheat;
		}

		if (!this->hardwireQ) {
			this->grid[i].Qimpur=QiSpline.get(log10(this->grid[i].P));
			if (this->grid[i].Qimpur < 0.0) this->grid[i].Qimpur=0.0;
		}

		if (this->output) 
			fprintf(fp, "%d %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg\n", i, this->grid[i].P, this->grid[i].rho, this->EOS->A[1]*(1.0-this->EOS->Yn), 
				this->EOS->Z[1], this->EOS->Yn,this->EOS->A[1],this->EOS->ptot(), Tmelt, GammaT/1e8, LoverT*1e8, this->grid[i].r, 
				1.0/sqrt(1.0-2.0*6.67e-8*2e33*this->mass/(9e20*this->grid[i].r)));
			//,AASpline.get(log10(this->grid[i].rho)),  ZZSpline.get(log10(this->grid[i].rho)), this->grid[i].Qimpur, this->grid[i].Qheat);
  	}

	if (this->output) fclose(fp);

	if (!this->hardwireQ) QiSpline.tidy();

	if (this->resume) read_T_profile_from_file();

  	printf("Grid has %d points, delx=%lg, Pb=%lg, rhob=%lg, Pt=%lg, rhot=%lg, thickness=%lg m\n", 
			this->N, this->dx, this->grid[this->N].P,this->grid[this->N].rho,this->grid[1].P,this->grid[1].rho,(this->grid[0].r-this->grid[this->N+1].r)*1e-2);
	if (!this->hardwireQ)
		printf("Total heat release is %lg MeV\n",Qtot);
}


void Crust::get_TbTeff_relation(void)
// reads in the Flux-T relation from the data file output by makegrid.cc
{
	double *temp, *flux;  // temporary storage to initialize the spline
	int npoints = 195;  //  needs to be >= number of points read in
	temp = new double [npoints+1];
	flux = new double [npoints+1];
	
	// the file "envelope_data/grid" is made by makegrid.cc
	// it contains  (column depth, T, flux)  in cgs
	FILE *fp;
	if (this->use_my_envelope) {
		if (this->EOS->B == 1e15) fp = fopen("envelope_data/grid_1e15_nopotek","r");
		else if (this->EOS->B == 1e14) fp = fopen("envelope_data/grid_1e14_potek","r");
		else if (this->EOS->B == 3e14) fp = fopen("envelope_data/grid_3e14_potek","r");
		else if (this->EOS->B == 3e15) fp = fopen("envelope_data/grid_3e15_potek","r");
		else {
			printf("Don't know which envelope model to use for this B!\n");
			exit(1);
		}
	} else {
		if (this->gpe) fp = fopen("envelope_data/grid_He4","r");
		else fp = fopen("envelope_data/grid_He9","r");
	}
	FILE *fp2=NULL;
	if (this->output) fp2=fopen("out/TbTeff", "w");
	
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
	
	delete [] temp;
	delete [] flux;
}


void Crust::set_temperature_profile(double *rhovec,double *Tvec,int nvec) 
// sets the temperature profile by interpolating between the specified (density, temperature) pairs
{
	// if density is <0 it means the base density
	// if density is 0 it means the top density
	// if temperature is <=0 it means the core temperature
	for (int i=1; i<nvec; i++) {
		if (Tvec[i] <= 0.0) Tvec[i]=this->Tc;
		if (rhovec[i] < 0.0) rhovec[i] = this->grid[this->N].rho;
		if (rhovec[i] == 0.0) {
			rhovec[i] = this->grid[1].rho;
		}
	}	
	if (rhovec[nvec-1] != this->grid[this->N].rho) {  // if we didn't specify it in the file,
									// set the temperature of the base to the core temperature
		nvec++;
		rhovec[nvec-1] = this->grid[this->N].rho;
		Tvec[nvec-1] = this->Tc;
	}

	// now assign initial temperatures to the grid
	// the temperature is linearly interpolated in log density
	for (int i=1; i<=this->N+1; i++) {
		double Ti;
		if (i==1) {
			Ti=Tvec[1]; 
			this->Tt=Ti;
		} else {
			if (i==this->N+1) {
				Ti=Tvec[nvec];		
			} else {
				int	j=0; 
				while (rhovec[j] < this->grid[i].rho && j<nvec) j++;
				Ti = pow(10.0,log10(Tvec[j-1]) + log10(Tvec[j]/Tvec[j-1])*log10(this->grid[i].rho/rhovec[j-1])/log10(rhovec[j]/rhovec[j-1]));
			}
		}	
		
		this->grid[i].T=Ti;
	}
}


void Crust::read_T_profile_from_file(void)
{
	int npoints;
	double dd, tt;

	FILE *fp = fopen("out/out","r");

	// first read the header
	fscanf(fp,"%d %lg\n",&npoints,&dd);
	if (npoints != this->N+1) {
		printf("Problem reading previous T profile: number of grid points is different\n");
		exit(0);
	}
	
	while (!feof(fp)) {	
		fscanf(fp,"%lg\n",&tt);		
		for (int i=1; i<=npoints; i++) {
			double T;
			fscanf(fp,"%lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg\n",&dd,&T,&dd,&dd,&dd,&dd,&dd,&dd,&dd,&dd,&dd,&dd,&dd);
			this->grid[i].T = T;
		}
	}

	fclose(fp);

	this->last_time_output = tt;	
	this->timesofar = tt / this->ZZ;
}







// --------------------------------- Time evolution ---------------------------------------------


void Crust::evolve(double timetorun, double mdot) {

	// time is the time to evolve in days
	// mdot is the accretion rate in Eddington units
	this->outburst_duration=timetorun/(365.0*this->ZZ);
	printf("Now evolve in time for %lg days at mdot=%lg (star time=%lg yrs)\n",timetorun,mdot,this->outburst_duration);
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
		this->ODE.set_bc(i,this->grid[i].T);
	}
	this->ODE.go(0.0, this->outburst_duration*3.15e7, this->outburst_duration*3.15e7*0.01,1e-7);
	stop_timing(&timer,"this->ODE.go");
	printf("number of steps = %d\n", this->ODE.kount);

	for (int i=1; i<=this->N+1; i++) {
		this->grid[i].T=ODE.get_y(i,this->ODE.kount);
	}

	// output results
	if (this->output) {
		printf("Starting output\n");
		if (this->last_time_output == 0.0) {
			this->fp=fopen("out/out","w");
		   	this->fp2=fopen("out/prof","w");
		} else {
			this->fp=fopen("out/out","a");
	   		this->fp2=fopen("out/prof","a");
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




void Crust::output_result_for_step(int j, FILE *fp, FILE *fp2,double timesofar,double *last_time_output) 
{
	// Output if enough time has elapsed
	if (fabs(log10(this->ODE.get_x(j)/this->ODE.get_x(j-1))) >= 0.01) {

		// get CP,K,eps,eps_nu at each point on the grid
		for (int i=1; i<=this->N+1; i++) {
			this->grid[i].T=this->ODE.get_y(i,j);
			calculate_vars(i);
		}
		// outer boundary
		double T0;
		this->grid[1].T=this->ODE.get_y(1,j);
		outer_boundary();
		T0=this->grid[0].T;

		// timestep
		double dt;
		if (j==1) dt=this->ODE.get_x(j); else dt=this->ODE.get_x(j)-this->ODE.get_x(j-1);

		// heat fluxes on the grid
		double *TT=new double [this->N+2];
		for (int i=1; i<=this->N+1; i++) TT[i]=this->ODE.get_y(i,j);
		double FF = calculate_heat_flux(1,TT);
		for (int i=1; i<=this->N+1; i++) this->grid[i].F = calculate_heat_flux(i,TT);
		delete [] TT;

		// total neutrino luminosity
		double Lnu=0.0;
		for (int i=1; i<=this->N; i++) Lnu += this->grid[i].NU*this->dx*this->grid[i].P/this->g;
	
		// we output time, fluxes and TEFF that are already redshifted into the observer frame
		// out/prof
		fprintf(fp2, "%lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg\n", (timesofar+this->ODE.get_x(j))*this->ZZ, 
			pow((this->radius/11.2),2.0)*this->grid[2].F/(this->ZZ*this->ZZ), pow((this->radius/11.2),2.0)*FF/(this->ZZ*this->ZZ),
			this->ODE.get_y(this->N-5,j), pow((this->g/2.28e14)*TEFF.get(this->ODE.get_y(1,j))/5.67e-5,0.25)/this->ZZ, 
			this->ODE.get_y(1,j), pow((this->g/2.28e14)*TEFF.get(this->ODE.get_y(1,j))/5.67e-5,0.25),
			pow((this->radius/11.2),2.0)*this->grid[this->N+1].F/(this->ZZ*this->ZZ),pow((this->radius/11.2),2.0)*this->grid[this->N].F/(this->ZZ*this->ZZ),
			4.0*M_PI*pow(1e5*this->radius,2.0)*Lnu/(this->ZZ*this->ZZ), dt);
			
		if ((fabs(log10(fabs(timesofar+this->ODE.get_x(j))*this->ZZ)-log10(fabs(*last_time_output))) >= 1000.0) ||
			(fabs(timesofar)+this->ODE.get_x(j))*this->ZZ < 1e10) {
			// temperature profile into out/out
			fprintf(fp,"%lg\n",this->ZZ*(timesofar+this->ODE.get_x(j)));
			for (int i=1; i<=this->N+1; i++)
				fprintf(fp, "%lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg\n", 
					this->grid[i].P, this->ODE.get_y(i,j), this->grid[i].F, this->grid[i].NU, this->g*(this->grid[i+1].F-this->grid[i].F)/(this->dx*this->grid[i].P), this->grid[i].rho, this->grid[i].CP*this->grid[i].rho, 
					this->ODE.get_d(i,j),1e8*pow(this->grid[i].P/2.521967e17,0.25), this->grid[i].K, 2.521967e-15*pow(this->ODE.get_y(i,j),4)/this->grid[i].P,
					this->grid[i].NU,this->grid[i].EPS);
 		}
		*last_time_output=(timesofar+this->ODE.get_x(j))*this->ZZ;

	}


}





// --------------------------------- Crust properties ---------------------------------------------


void Crust::set_composition(void)
// sets the composition appropriate for the current pressure
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

	this->CP_grid = matrix(this->N+2,this->nbeta);	
	this->K1_grid = matrix(this->N+2,this->nbeta);	
	this->K0_grid = matrix(this->N+2,this->nbeta);	
	this->KAPPA_grid = matrix(this->N+2,this->nbeta);	
	this->K1perp_grid = matrix(this->N+2,this->nbeta);	
	this->K0perp_grid = matrix(this->N+2,this->nbeta);	
	this->NU_grid = matrix(this->N+2,this->nbeta);	
	this->EPS_grid = matrix(this->N+2,this->nbeta);	

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
	if (EOS->B > 0.0) sprintf(s,"out/precalc_results_%lg",log10(EOS->B));
	else sprintf(s,"out/precalc_results_0");
	if (!this->force_precalc) fp=fopen(s,"r");
	// if unsuccessful (or if precalc is set) we need to recalculate
	if (fp == NULL) {
		if (this->output)
			fp=fopen(s,"w");
		printf("Precalculating quantities and writing to file %s...\n",s);

		for (int i=1; i<=this->N+1; i++) {
	
			EOS->P=this->grid[i].P;
			EOS->rho = this->grid[i].rho;
			set_composition();
			
			if (this->output) 
				fprintf(fp, "Grid point %d  P=%lg  rho=%lg  A=%lg  Z=%lg Yn=%lg:  T8,CP,K,eps_nu,eps_nuc\n",
					i, this->grid[i].P, this->grid[i].rho, (1.0-EOS->Yn)*EOS->A[1], EOS->Z[1], EOS->Yn);
		
			double heating = crust_heating(i);
		
			for (int j=1; j<=this->nbeta; j++) {		
				double beta = this->betamin + (j-1)*(this->betamax-this->betamin)/(1.0*(this->nbeta-1));
				EOS->T8 = 1e-8*pow(10.0,beta);

				if (i == this->N+1) {
					this->CP_grid[i][j] = this->C_core * EOS->T8;
					this->NU_grid[i][j] = this->Lnu_core_norm * pow(EOS->T8, Lnu_core_alpha);
					this->EPS_grid[i][j] = 0.0; // no core heating
					this->K0_grid[i][j]=this->K0_grid[i-1][j];
					this->K1_grid[i][j]=this->K1_grid[i-1][j];
					this->K1perp_grid[i][j]=this->K1perp_grid[i-1][j];
					this->K0perp_grid[i][j]=this->K0perp_grid[i-1][j];
					
				} else {
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
					this->K0_grid[i][j]=EOS->rho*Kcond/this->grid[i].P;
					this->K0perp_grid[i][j]=EOS->rho*Kcondperp/this->grid[i].P;

					EOS->Qimp=1.0;
					//Kcond = EOS->K_cond(EOS->Chabrier_EF());
					//Kcondperp=Kcond;
					Kcond = EOS->potek_cond();
					Kcondperp = EOS->Kperp;
					this->K1_grid[i][j]=EOS->rho*Kcond/this->grid[i].P;
					this->K1perp_grid[i][j]=EOS->rho*Kcondperp/this->grid[i].P;

					EOS->Qimp=Q_store;  // restore to previous value

					// conductivity due to radiation
					(void) EOS->opac();  // call to kappa sets the variable kappa_rad
					this->KAPPA_grid[i][j] = 3.03e20*pow(EOS->T8,3)/(EOS->kappa_rad*this->grid[i].P);

				}

				if (this->output) 
					fprintf(fp, "%lg %lg %lg %lg %lg %lg %lg %lg %lg\n", EOS->T8, this->CP_grid[i][j], 
						this->K0_grid[i][j],this->K1_grid[i][j], this->K0perp_grid[i][j],this->K1perp_grid[i][j],
						this->NU_grid[i][j], this->EPS_grid[i][j], this->KAPPA_grid[i][j] );
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
			}
		}
		fclose(fp);
			
	}
	
}

double Crust::crust_heating(int i) 
// calculates the crust heating for grid point i
// units are erg/g/s  divided by (mdot*g)
// (the mdot*g factor is put back in when we calculate dTdt, so that we don't need to precalculate when changing either mdot or g)
{
	double eps=0.0,P = this->grid[i].P;

	// if we are heating on < 1 day timescale then it is a magnetar
	if (this->outburst_duration<1.0/365.0) {
		
		// eps in erg/g/s
		double eps_heat = 1e25/(this->grid[i].rho*this->outburst_duration*3.15e7);
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
		
		{ // the above assumed 1e25 erg/cm^3 deposited energy; now apply a multiplier as specified in the inlist.dat
			double ener;
			if (this->grid[i].rho>4e11) ener = this->energy_deposited_inner;
			else ener = this->energy_deposited_outer;
			ener *= pow(this->grid[i].rho/1e10,this->energy_slope);
			eps *= ener;
		}
		
	} else {  // otherwise we are doing an accreting neutron star

		if (!this->hardwireQ) {   // the profile of Q(rho) was specified in the crust model
			eps = this->grid[i].Qheat*8.8e4*9.64e17/(this->grid[i].P*this->dx);
		} else {
			// simple "smeared out" heating function, 1.5MeV in inner crust, 0.2MeV in outer crust (as in BC09)
			eps += eps_from_heat_source(P,1e16,1e17,1.5);	
			eps += eps_from_heat_source(P,3e12,3e15,0.2);

			// Extra heat source in the ocean
			if (this->extra_heating) {	
				// Put all of the extra heat into one grid point
				//if (this->grid[i].P*exp(-0.5*this->dx) <this->extra_y*2.28e14 && this->grid[i].P*exp(0.5*this->dx)>this->extra_y*2.28e14)
				//		eps+=8.8e4*this->extra_Q*9.64e17/(P*this->dx);
				double heating_spread = 3.0;
				eps += eps_from_heat_source(P,this->extra_y/heating_spread,this->extra_y*heating_spread,this->extra_Q);				
			}
		}
	}

	return eps;	
}



double Crust::eps_from_heat_source(double P,double y1,double y2,double Q_heat)
// returns the local heating rate at pressure P for a heat source distributed between effective pressures y1 and y2
{
	double eps = 0.0;
	double P1 = P*exp(-0.5*this->dx);
	double P2 = P*exp(0.5*this->dx);
	double geff=2.28e14;

	if (P1 > y1*geff && P2 < y2*geff)   // we are within the heating zone
		eps=8.8e4*Q_heat*9.64e17/(P*log(y2/y1));
	if (P1 < y1*geff && P2 < y2*geff && y1*geff<P2) {   // left hand edge of heated region
		eps=8.8e4*Q_heat*9.64e17/(P*log(y2/y1));
		eps *= log(P2/(y1*geff))/this->dx;
	}
	if (P1 > y1*geff && P2 > y2*geff && y2*geff>P1) {  // right hand edge of heated region
		eps=8.8e4*Q_heat*9.64e17/(P*log(y2/y1));
		eps *= log(y2*geff/P1)/this->dx;
	}
	
	return eps;
}




// ----------------------------- Integration -------------------------------------------------------------


void Crust::derivs(double t, double T[], double dTdt[])
// calculates the time derivatives for the whole grid
{
	// First calculate quantities at each grid point
	for (int j=1; j<=this->N+1; j++) {
		this->grid[j].T=T[j];
		calculate_vars(j);
	}
	outer_boundary();
	T[0]=this->grid[0].T;

	// determine the fluxes at the half-grid points
	//  this->grid[i].F is the flux at i-1/2
  	for (int i=1; i<=this->N+1; i++)   this->grid[i].F = calculate_heat_flux(i,T);	
	
	// Calculate the derivatives dT/dt
	for (int i=1; i<=this->N; i++) {
  		dTdt[i]=this->g*pow(this->grid[0].r/this->grid[i].r,4.0)*(this->grid[i+1].F-this->grid[i].F)/(this->dx*this->grid[i].CP*this->grid[i].P);
		if (this->nuflag) dTdt[i]+=-(this->grid[i].NU/this->grid[i].CP);
		if (this->accreting) dTdt[i]+=this->grid[i].EPS/this->grid[i].CP;
	}
	// the cell at N+1 represents the core
  	dTdt[this->N+1] = (-this->grid[this->N+1].F * 4.0*M_PI*pow(1e5*this->radius,2.0) - this->grid[this->N+1].NU) / this->grid[this->N+1].CP;
}

double Crust::calculate_heat_flux(int i, double *T)
{
	double flux;
	if (i>1 || (this->accreting && this->outburst_duration > 1.0/365.0 && !this->force_cooling_bc))
//		if (i>1 || (this->accreting && EOS->B == 0.0))   
		// use this inside the grid, or at the surface when we are accreting (which 
		// fixes the outer temperature)
		if (i==this->N+1)
			flux = this->grid[i-1].K*(T[i]-T[i-1])/this->dx;	
		else
			flux = 0.5*(this->grid[i].K+this->grid[i-1].K)*(T[i]-T[i-1])/this->dx;	
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
	for (int j=k; j<=k2; j++) {
		this->grid[j].T=T[j];
		calculate_vars(j);
	}
	if (i==1) {
		this->grid[1].T=T[1];
		outer_boundary();
		T[0]=this->grid[0].T;
	}

	double f;
	if (i<this->N+1) {
		f=this->g*pow(this->grid[0].r/this->grid[i].r,4.0)*(calculate_heat_flux(i+1,T)-calculate_heat_flux(i,T))/(this->dx*this->grid[i].CP*this->grid[i].P);
		if (this->nuflag) f+=-(this->grid[i].NU/this->grid[i].CP);	
		if (this->accreting) f+=this->grid[i].EPS/this->grid[i].CP;
	} else {
		f = (-calculate_heat_flux(i,T) * 4.0*M_PI*pow(1e5*this->radius,2.0) - this->grid[i].NU) / this->grid[i].CP;
	}
	
	return f;
}

void Crust::outer_boundary(void)
{
	if (this->accreting && this->Tt>0.0) this->grid[0].T=this->Tt;   // constant temperature during accretion
	else this->grid[0].T=this->grid[1].T*(8.0-this->dx)/(8.0+this->dx);   // assumes radiative zero solution, F\propto T^4
	this->grid[0].K=this->grid[1].K; this->grid[0].CP=this->grid[1].CP;
	if (this->nuflag) this->grid[0].NU=this->grid[1].NU; else this->grid[0].NU=0.0;
	if (this->accreting) this->grid[0].EPS=this->grid[1].EPS; else this->grid[0].EPS=0.0;
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

  {
	int i=n;
	T[i]*=1.0+e; f=dTdt(i,T);
	T[i]/=1.0+e; dfdT[i][i]=(f-dfdt[i])/(T[i]*e);
	T[i-1]*=1.0+e; f=dTdt(i,T);
	T[i-1]/=1.0+e; dfdT[i][i-1]=(f-dfdt[i])/(T[i-1]*e);
  }

}  


void Crust::calculate_vars(int i)
{
	double T=this->grid[i].T;
	double P=this->grid[i].P;
	double *CP=&this->grid[i].CP;
	double *K=&this->grid[i].K;
	double *NU=&this->grid[i].NU;
	double *EPS=&this->grid[i].EPS;

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
		if (this->grid[i].rho > this->Qrho) Qval=this->Qinner; else Qval=EOS->Qimp;
//		if (P>2.28e29) Qval=this->Qinner; else Qval=EOS->Qimp;
	} else {
		Qval = this->grid[i].Qimpur;	
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





