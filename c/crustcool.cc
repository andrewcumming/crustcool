// crustcool.cc
//
// A time dependent code to model the thermal evolution
// of an accreting neutron star crust
// based on "cool.cc" which makes long X-ray burst lightcurves

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <time.h>
#include "../h/nr.h"
#include "../h/nrutil.h"
#include "../h/odeint.h"
#include "../h/eos.h"
#include "../h/spline.h"

#pragma mark ====== Declarations ======
double dTdt(int i, double *T);
void derivs(double t, double T[], double dTdt[]);
void jacobn(double, double *, double *, double **, int);
void calculate_vars(int i, double T, double y, double *CP, double *K, double *NU, double *EPS);
void outer_boundary(double T1, double K1, double CP1, double NU1, double EPS1, double *T0, double *K0, double *CP0, double *NU0, double *EPS0);
void inner_boundary(double TN, double KN, double CPN, double NUN, double EPSN, double *TN1, double *KN1, double *CPN1, double *NUN1, double *EPSN1);
double calculate_heat_flux(int i, double *T);
void set_up_initial_temperature_profile(void);
void set_up_initial_temperature_profile_piecewise(char *fname);
void precalculate_vars(void);
void set_up_grid(int ngrid,const char*fname);
void get_TbTeff_relation(void);
double crust_heating_rate(double y);
double crust_heating(int i);
double energy_deposited(int i);
void output_result_for_step(int j, FILE *fp, FILE *fp2, double timesofar);
void start_timing(clock_t *timer);
void stop_timing(clock_t *timer, const char*string);
void read_in_data(const char *fname);
void calculate_cooling_curve(void);
double calculate_chisq(void);
void set_composition(void);
void set_ns_parameters(double mass, double radius);
void set_ns_redshift(double R);
extern "C"{
  void condegin_(double *temp,double *densi,double *B,double *Zion,double *CMI,double *CMI1,double *Zimp, double *RSIGMA,double *RTSIGMA,double *RHSIGMA,double *RKAPPA,double *RTKAPPA,double *RHKAPPA);
   }
double potek_cond(double *Kperp);
void output_some_EOS_info();
void heatderivs(double t, double T[], double dTdt[]);
void heatderivs2(double t, double T[], double dTdt[]);
void make_TbTeff_relation(void);
void surfderivs(double y, double T[], double dTdy[]);
double find_surf_eqn(double y);

// global variables
struct globals {
  	int N;  
  	double dx; 
  	double *P, *CP, *K, *F, *NU, *rho, *EPS, *Tmelt, *GammaT, *LoverT;
	double *Qheat,*Qimp;
	double **CP_grid, **K1_grid, **K0_grid, **NU_grid, **EPS_grid, **K1perp_grid, **K0perp_grid;
	double betamin, betamax, deltabeta;
	int nbeta;
  	double g, ZZ,mass, radius;
  	int nuflag;  
	double time_to_run;
	double Pb, Pt, yt;
	double rhot,rhob;
	double mdot;
	double Tt, Fin, Tc;
	int accreting;
	double *obs_time, *obs_temp, *obs_err;
	int obs_n;
	int output;
	int include_latent_heat;
	int running_from_command_line;
	int hardwireQ;
	double Qrho;
	int use_piecewise;
	int force_precalc;
	int include_sph;
	FILE *fp,*fp2;
	double Qinner;
	double energy_deposited_outer;
	double energy_deposited_inner;
	double outburst_duration;
	int instant_heat;
	double surfF, surfy;
	double angle_mu;
	int gpe;
} G;

Ode_Int ODE;
Eos EOS;
Spline RHO;
Spline TEFF;
Spline AASpline; 
Spline ZZSpline;
Spline YnSpline;


#pragma mark =========== Code ============

int main(int argc, char *argv[])
{
	double F, y, K, CP, Kh, Th, Ttop, Fb, mass;
	clock_t timer;
	
	int do_mcmc=0;
	
	long idum=-32094034;

	// initialize EOS routines
  	EOS.init(1); 
  	EOS.X[1] = 1.0;   // we only have one species at each depth;

	//output_some_EOS_info();  

	// ----------------------------------------------------------------------
 	//   Set parameters

	// first, some defaults	
	mass=1.62;
	G.radius=11.2;
	int ngrid=100;
	G.include_latent_heat=0;
	G.nuflag = 1;
	G.force_precalc=0;
	G.Qinner=-1.0;
	G.running_from_command_line=1;	
	G.outburst_duration = (1.0/24.0) * 1.0/(365.0);  // rapid heating for magnetar case	
	EOS.accr = 0;   // set crust composition
	EOS.B=0; //2.2e14;   // magnetic field in the crust   (set B>0 for magnetar case)
	EOS.gap = 1;    // 0 = no gap, normal neutrons
	EOS.kncrit = 0.0;  // neutrons are normal for kn<kncrit (to use this set gap=4)
	G.mdot=0.1;
	G.energy_deposited_outer=1.0;
	G.energy_deposited_inner=-1.0;
	do_mcmc=0;	
	G.rhob=1e15; 
	G.rhot=1e3;
	G.use_piecewise=0;
	G.Qrho=1e12;
	G.instant_heat = 0;
	G.include_sph=1;
	G.angle_mu=-1.0;
	G.gpe=0;
	G.yt=1e12;
	
	// now read from the file 'init.dat'
	char fname[40];
	char fnamedefault[10]="init.dat";
	if (argc >1) {
		strcat(fname,"init/init.dat.");
		strcat(fname,argv[1]);
	} else {
		strcat(fname,fnamedefault);
	}
	printf("============================================\n");
	printf("Reading input data from %s\n",fname);
	FILE *fp = fopen(fname,"r");
	char s1[100];
	char s[100];
	double x;				
	while (!feof(fp)) {   // we read the file line by line
		char *e=fgets(s1,200,fp);		
		// ignoring lines that begin with \n (blank) or with # (comments)
		// or with $ (temperature profile)
		if (strncmp(s1,"#",1) && strncmp(s1,"\n",1) && strncmp(s1,">",1)) {
			sscanf(s1,"%s\t%lg\n",s,&x);
			if (!strncmp(s,"Bfield",6)) EOS.B=x;
			if (!strncmp(s,"Tc",2)) G.Tc=x;
			if (!strncmp(s,"Tt",2)) G.Tt=x;
			if (!strncmp(s,"SFgap",5)) EOS.gap=(int) x;
			if (!strncmp(s,"ngrid",5)) ngrid=(int) x;
			if (!strncmp(s,"kncrit",6)) EOS.kncrit=x;
			if (!strncmp(s,"mdot",4)) G.mdot=x;
			if (!strncmp(s,"mass",4)) mass=x;
			if (!strncmp(s,"sph",3)) G.include_sph=(int) x;
			if (!strncmp(s,"gpe",3)) G.gpe=(int) x;
			if (!strncmp(s,"radius",6)) G.radius=x;
			if (!strncmp(s,"Edep",4)) G.energy_deposited_outer=x;
			if (!strncmp(s,"ytop",4)) G.yt=x;
			if (!strncmp(s,"Einner",6)) G.energy_deposited_inner=x;
			if (!strncmp(s,"Qimp",4)) EOS.Q=x;
			if (!strncmp(s,"Qrho",4)) G.Qrho=x;
			if (!strncmp(s,"rhob",4)) G.rhob=x;
			if (!strncmp(s,"rhot",4)) G.rhot=x;
			if (!strncmp(s,"precalc",7)) G.force_precalc=(int) x;
			if (!strncmp(s,"instant",7)) G.instant_heat=(int) x;
			if (!strncmp(s,"Qinner",6)) G.Qinner=x;
			if (!strncmp(s,"timetorun",9)) G.time_to_run=24.0*3600.0*x;			
			if (!strncmp(s,"toutburst",9)) G.outburst_duration=x;
			if (!strncmp(s,"piecewise",9)) G.use_piecewise=(int) x;
			if (!strncmp(s,"neutrinos",9)) G.nuflag=(int) x;
			if (!strncmp(s,"accreted",8)) EOS.accr=(int) x;
			if (!strncmp(s,"angle_mu",8)) G.angle_mu=x;
			if (!strncmp(s,"latent_heat",11)) G.include_latent_heat=(int) x;
		}
	}

	fclose(fp);

	if (G.Qinner == -1.0) G.Qinner=EOS.Q;
	if (G.energy_deposited_inner == -1.0) G.energy_deposited_inner = G.energy_deposited_outer;
		
	//	read_in_data("data/1731");  // READ IN observed lightcurve
		read_in_data("data/1659");  // READ IN observed lightcurve
	//	read_in_data("data/XTEJ");  // READ IN observed lightcurve
	//	read_in_data("data/terz");  // READ IN observed lightcurve
	//	read_in_data("data/terz2");  // READ IN observed lightcurve
	//	read_in_data("data/0748");  // READ IN observed lightcurve

	set_ns_parameters(mass,G.radius);
	G.outburst_duration /= G.ZZ;   // redshift the outburst duration (shorter time on the NS surface)

 	// set up the hydrostatic grid
	set_up_grid(ngrid,"data/crust_model_shell");

	// ----------------------------------------------------------------------------------

	//make_TbTeff_relation(); 

	// read in Tb-Teff relation
	get_TbTeff_relation();
	
  	// precalculate CP, K, eps_nu as a function of T on the grid
	// note that for K this is done in such a way that we do not need to recalculate
	// if Qimp changes
	start_timing(&timer);
  	precalculate_vars();
	stop_timing(&timer,"precalculate_vars");

	
  	// initialize the integrator
  	ODE.init(G.N+1);
  	ODE.verbose=0;   // print out each timestep if this is set to 1
  	ODE.stiff=1; ODE.tri=1;  // stiff integrator with tridiagonal solver
  
	// do we want output or not?
	G.output=!do_mcmc;

 	// set up the initial temperature profile
  	if (G.use_piecewise) set_up_initial_temperature_profile_piecewise(fname);
	else set_up_initial_temperature_profile();

	if (G.output) {
  		G.fp=fopen("gon_out/out","w");
  		G.fp2=fopen("gon_out/prof","w");
  		fprintf(G.fp,"%d\n",G.N+1);
	}

	// calculate the cooling curve
	calculate_cooling_curve();

	// get chisq
	double chisq = calculate_chisq();


	if (do_mcmc) {
		
		//G.output=0;  // switch off output to other files (for speed)
  		FILE *fp_out=fopen("out/settle2.mcmc","w");	// output MCMC chain
		double Q, prob;
		int accept=0;
		int count=0, accept_count=0;
		
		if (EOS.Q<1e-3 || G.Qinner<1e-3) {
			printf("One or both Q values are too small: setting to 1e-3\n");
			if (EOS.Q<1e-3) EOS.Q=1e-3;
			if (G.Qinner < 1e-3) G.Qinner=1e-3;
		}
		
		double Tc_old=G.Tc;
		double Tt_old=G.Tt;
		double Q_old=log10(EOS.Q);
		double Qinner_old=log10(G.Qinner);
		double M_old = mass;
		double R_old=G.radius;
		double chisq_old = chisq;	
		double mdot_old = G.mdot;
		
		while (1) {

			double f=1.0;
		 	// take a step 
		while (Q=Q_old + f*0.1*gasdev(&idum), Q<-3.0) {};
//		while (Q=Q_old + f*0.02*gasdev(&idum), Q<-3.0) {};
//			while (G.Qinner=Qinner_old + 0.2*gasdev(&idum), G.Qinner<-3.0) {};

//   		while (G.Tc=Tc_old + f*3e6*gasdev(&idum), G.Tc < 2e7);

   		while (G.Tt=Tt_old + f*3e7*gasdev(&idum), G.Tt < 0.0);

//   		while (G.mdot=mdot_old + f*0.03*gasdev(&idum), G.mdot < 0.0 || G.mdot > 0.5);

   		//while (mass=M_old + f*0.1*gasdev(&idum), mass < 1.0 || mass >3.0);
   		//while (radius=R_old + f*0.3*gasdev(&idum), radius < 9.0 || radius > 16.0);
			EOS.Q=pow(10.0,Q);
			//G.Qinner=pow(10.0,G.Qinner);
			G.Qinner=EOS.Q;

			// recalculate the model
			set_ns_parameters(mass,G.radius);
			//set_up_grid(ngrid,"data/crust_model");
			//precalculate_vars();
			set_up_initial_temperature_profile();
			calculate_cooling_curve();
			chisq=calculate_chisq();
			prob=exp((chisq_old-chisq)/2.0);

   		// adjust probabilities to account for the prior
			//   prob*=Q/Qnew;
   		prob*=Tt_old/G.Tt;
   		prob*=Tc_old/G.Tc;

			count++;
			accept=0;
			if (prob>=1.0) accept=1; else if(ran2(&idum) < prob) accept=1;
			if (accept) {
				accept_count++;
   			fprintf(fp_out, "%lg %lg %lg %lg %lg %lg %lg %lg %lg %lg\n", 
						G.Tc, G.Tt, EOS.Q, mass, G.radius, G.ZZ, G.g*1e-14, chisq, G.mdot, G.Qinner);
				fflush(fp_out);
				Tc_old=G.Tc; Tt_old=G.Tt; Q_old=log10(EOS.Q); chisq_old = chisq; 
				M_old = mass; R_old = G.radius; mdot_old = G.mdot; Qinner_old = log10(G.Qinner);
				printf("accepted:");
			} else {
				printf("not accepted:");
			}
			printf("%lg %lg %lg %lg %lg %lg %lg %lg %lg\n---------------------------\n",
			 		G.Tc, G.Tt, EOS.Q, G.Qinner, mass, G.radius, G.mdot, chisq, 1.0*accept_count/count );
		}
	}

/*
	if (do_mcmc==2) {
		
		double chisq_old, Tc_old, Tt_old, Q_old, Q, prob;
		int accept=0;
		int count=0, accept_count=0;

  		FILE *fp_out=fopen("out/settle2.mcmc","w");

		double *oldT;
		if (1) {
			fprintf(fp_out,"%d\n",G.N);
			oldT=vector(1,G.N);
			for (int i=1; i<=G.N; i++) {
				fprintf(fp_out, "%lg ", G.y[i]);
			}
			fprintf(fp_out,"\n");
			for (int i=1; i<=G.N; i++) {
				oldT[i]=ODE.get_y(i,1);
				fprintf(fp_out, "%lg ", oldT[i]);
			}
			fprintf(fp_out,"\n");
		}

		double grad_term=0.0;
		for (int i=2; i<=G.N; i++) {
			grad_term += pow((oldT[i]-oldT[i-1])/(oldT[i]+oldT[i-1]),2.0);
		}
		chisq*=(1.0+grad_term);
		Tc_old=G.Tc; Tt_old=G.Tt; Q_old=log10(EOS.Q); chisq_old = chisq;

		
		
		while (1) {

		 	// take a step and adjust the initial profile
			if (0) {
		
				while (Q=Q_old + 0.07*gasdev(&idum), Q<-3.0) {};
   			while (G.Tc=Tc_old + 0.7e6*gasdev(&idum), G.Tc < 0.0);
   			while (G.Tt=Tt_old + 0.7e7*gasdev(&idum), G.Tt < 0.0);
				EOS.Q=pow(10.0,Q);
				set_up_initial_temperature_profile();
			
			} else {
	
				for (int i=1; i<=G.N; i++) {
					double newT;
					while (newT=oldT[i]*(1.0+0.12*gasdev(&idum)), newT<0.2*oldT[i]) {};
					ODE.set_bc(i,newT);
				}
			}		
		
			calculate_cooling_curve();
			chisq=calculate_chisq();
			grad_term=0.0;
			for (int i=2; i<=G.N; i++) {
				grad_term += pow((oldT[i]-oldT[i-1])/(oldT[i]+oldT[i-1]),2.0);
			}
			chisq*=(1.0+grad_term);
			prob=exp((chisq_old-chisq)/2.0);
   		// adjust probabilities to account for the prior
			//   prob*=Q/Qnew;
   		prob*=Tt_old/G.Tt;
   		prob*=Tc_old/G.Tc;

			count++;
			accept=0;
			if (prob>=1.0) accept=1; 
			else 
				if(ran2(&idum) < prob) accept=1;
			if (accept) {
				accept_count++;
   			//fprintf(fp_out, "%lg %lg %lg %lg %lg %lg %lg %lg %lg\n", G.Tc, G.Tt, EOS.Q, 10.0, 1.4, 1.31, G.g*1e-14, chisq, G.mdot);
				//fflush(fp_out);
				//Tc_old=G.Tc; Tt_old=G.Tt; Q_old=log10(EOS.Q); chisq_old = chisq;
				for (int i=1; i<=G.N; i++) {
					oldT[i]=ODE.get_y(i,1);
					fprintf(fp_out, "%lg ", oldT[i]);
				}
				fprintf(fp_out,"\n");
				fflush(fp_out);
				printf("accept: %lg %lg %lg %lg %lg\n", G.Tc, G.Tt, EOS.Q, chisq, 1.0*accept_count/count );
			} else {
				printf("not accepted: %lg %lg %lg %lg %lg\n", G.Tc, G.Tt, EOS.Q, chisq, 1.0*accept_count/count );
			}
			printf("---------------------------\n");
		}
	}

*/

	// tidy up
	if (G.output) {
		fclose(G.fp);
  		fclose(G.fp2);
	}
  	ODE.tidy(); 
	EOS.tidy();
	RHO.tidy();
  	free_vector(G.rho,0,G.N+2);
  	free_vector(G.CP,0,G.N+2);
  	free_vector(G.P,0,G.N+2);
  	free_vector(G.K,0,G.N+2);
  	free_vector(G.F,0,G.N+2);
  	free_vector(G.NU,0,G.N+2);
  	free_vector(G.EPS,0,G.N+2);
}



void calculate_cooling_curve(void) 
{
  	int nsteps=0;
  	double timesofar=0.0;
  	double dtstart=1e-6;  // start timestep in sec

  	while (true) {	
    	if (!G.running_from_command_line) {
			printf("Time to evolve (s) (enter 0 to end and write out summary info)? "); 
			scanf("%lg", &G.time_to_run);
		} else {
			if (timesofar>0.0) break;
			printf("Running for time %lg seconds\n", G.time_to_run);
		}
    	if (G.time_to_run==0.0) break;
		if (timesofar > 0.0) {   
			// set boundary condition for this integration using the last integration
			for (int i=1; i<=G.N+1; i++) ODE.set_bc(i, ODE.get_y(i,ODE.kount));
			dtstart=1e-6*G.time_to_run;
			if (EOS.B > 0.0) dtstart=1e3;
		}	
		G.accreting=0;
    
		// call the integrator
		clock_t timer;
		start_timing(&timer);
		ODE.dxsav=1e3;   // so that we get output starting at early enough times
		//ODE.go_scale(0.0, G.time_to_run, 1e4, derivs);
		ODE.go(0.0, G.time_to_run, dtstart, 1e-6, derivs);
		stop_timing(&timer,"ODE.go");
		
		// output 
		if (G.output) {
			for (int j=1; j<=ODE.kount; j++) output_result_for_step(j,G.fp,G.fp2,timesofar);
			fflush(G.fp); fflush(G.fp2);
		}
      timesofar+=G.time_to_run;
      nsteps+=ODE.kount;
      printf("number of steps = %d (total=%d) (time=%lg)\n", ODE.kount, nsteps, timesofar);

	}
}

double calculate_chisq(void)	
// uses the result of the cooling to calculate chi-squared
{
	// set up a spline which has the prediction for observed Teff vs time 
	Spline TE;
	double *yy = vector(1,ODE.kount);
	double *xx = vector(1,ODE.kount);
	for (int k=1; k<=ODE.kount; k++) { 
		xx[k]=ODE.get_x(k)*G.ZZ/(3600.0*24.0);
		yy[k]=1.38e-16*pow((TEFF.get(ODE.get_y(1,k))*(G.g/2.28e14))/5.67e-5,0.25)/(1.6e-12*G.ZZ);
	}
	TE.minit(xx,yy,ODE.kount);
	free_vector(xx,1,ODE.kount);
	free_vector(yy,1,ODE.kount);

	// calculate chisq
	double chisq=0.0;
	for (int i=1; i<=G.obs_n; i++) {
		chisq += pow((G.obs_temp[i] - TE.get(G.obs_time[i]))/G.obs_err[i],2.0);	
		//printf("%lg %lg %lg\n", G.obs_time[i], G.obs_temp[i], TE.get(G.obs_time[i]));
	}
	TE.tidy();
	printf("chisq = %lg/(%d-3) = %lg\n", chisq, G.obs_n, chisq/(G.obs_n-3));
	return chisq;

}



void read_in_data(const char *fname) 
{
	
	printf("Reading data from %s\n", fname);
	
	FILE *fp = fopen(fname,"r");
	
	double t0;
	fscanf(fp, "%lg %d\n", &t0, &G.obs_n);
	
	G.obs_time = vector(1,G.obs_n);
	G.obs_temp = vector(1,G.obs_n);
	G.obs_err = vector(1,G.obs_n);
		
	for (int i=1; i<=G.obs_n; i++) {
		double dummy;
		fscanf(fp,"%lg %lg %lg %lg %lg\n",  &G.obs_time[i], &dummy, &dummy, &G.obs_temp[i], 
					&G.obs_err[i]);
		G.obs_time[i]-=t0;
	}
	
	fclose(fp);	
}



void start_timing(clock_t *time)
{
	*time=clock();
}

void stop_timing(clock_t *time, const char*string)
{
  	printf("Time taken for %s =%lg s\n", string, (double) (clock()-*time)/((double) CLOCKS_PER_SEC)); 	
}


void output_result_for_step(int j, FILE *fp, FILE *fp2,double timesofar) 
{
	for (int i=1; i<=G.N+1; i++) calculate_vars(i,ODE.get_y(i,j),G.P[i],&G.CP[i],&G.K[i],&G.NU[i],&G.EPS[i]);
	double T0;
	outer_boundary(ODE.get_y(1,j),G.K[1],G.CP[1],G.NU[1],G.EPS[1],&T0,&G.K[0],&G.CP[0],&G.NU[0],&G.EPS[0]);
	for (int i=2; i<=G.N+1; i++) G.F[i]=0.5*G.g*(G.K[i]+G.K[i-1])*(ODE.get_y(i,j)-ODE.get_y(i-1,j))/G.dx;
	G.F[1]=0.5*(G.K[0]+G.K[1])*(ODE.get_y(1,j)-T0)/G.dx;
	double dt;
	if (j==1) dt=ODE.get_x(j); else dt=ODE.get_x(j)-ODE.get_x(j-1);


	double *TT;
	TT=vector(1,G.N+1);
	for (int i=1; i<=G.N+1; i++) TT[i]=ODE.get_y(i,j);
	double FF = calculate_heat_flux(1,TT);

	//FF=G.F[1];

	G.F[2] = calculate_heat_flux(2,TT);

	free_vector(TT,1,G.N+1);


	// output time, fluxes and TEFF that are already redshifted
	fprintf(fp2, "%lg %lg %lg %lg %lg %lg %lg\n", (timesofar+ODE.get_x(j))*G.ZZ, 
			G.F[2]/(G.ZZ*G.ZZ), FF/(G.ZZ*G.ZZ),
			ODE.get_y(G.N-5,j), pow((G.g/2.28e14)*TEFF.get(ODE.get_y(1,j))/5.67e-5,0.25)/G.ZZ, 
			ODE.get_y(1,j), pow((G.g/2.28e14)*TEFF.get(ODE.get_y(1,j))/5.67e-5,0.25));

	if (j % 1 == 0 || j==ODE.kount) {   // output every nth cycle
		// temperature profile
		fprintf(fp,"%lg\n",G.ZZ*(timesofar+ODE.get_x(j)));
		double del;
		for (int i=1; i<=G.N+1; i++) {      
			EOS.P=G.P[i]; EOS.T8=1e-8*ODE.get_y(i,j); EOS.rho=G.rho[i];//EOS.find_rho();
			if (i>1) del = (ODE.get_y(i+1,j)-ODE.get_y(i-1,j))/(2.0*G.dx*ODE.get_y(i,j));
			else del=1.0;
			fprintf(fp, "%lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg\n", 
				G.P[i], ODE.get_y(i,j), G.F[i], G.NU[i],
				G.g*(G.F[i+1]-G.F[i])/(G.dx*G.P[i]), EOS.rho, EOS.CP()*EOS.rho, 
				ODE.get_d(i,j),1e8*pow(EOS.P/2.521967e17,0.25), EOS.opac(), del, EOS.del_ad(),
				2.521967e-15*pow(ODE.get_y(i,j),4)/G.P[i]);	
 		}
	}
}




// --------------------------------------------------------------------------
// derivatives for the time-dependent code

void derivs(double t, double T[], double dTdt[])
// calculates the time derivatives for the whole grid
{
	//	printf("&& %lg\n", ODE.timestep);
  for (int j=1; j<=G.N; j++) calculate_vars(j,T[j], G.P[j], &G.CP[j], &G.K[j], &G.NU[j],&G.EPS[j]);
  outer_boundary(T[1],G.K[1],G.CP[1],G.NU[1],G.EPS[1],&T[0],&G.K[0],&G.CP[0],&G.NU[0],&G.EPS[0]);
  inner_boundary(T[G.N],G.K[G.N],G.CP[G.N],G.NU[G.N],G.EPS[G.N],
				&T[G.N+1],&G.K[G.N+1],&G.CP[G.N+1],&G.NU[G.N+1],&G.EPS[G.N+1]);
  for (int i=1; i<=G.N+1; i++)   G.F[i] = calculate_heat_flux(i,T);	

	if (G.include_latent_heat) {
  		for (int i=2; i<=G.N; i++) {
			double Gm = G.GammaT[i-1]/(0.5*(T[i-1]+T[i]));
			double Gp = G.GammaT[i]/(0.5*(T[i]+T[i+1]));
			if (Gm < 175.0 && Gp > 175.0) {
				double dxdT = (4.0/T[i]) * pow(175.0/Gm,4.0)/ (pow(Gp/Gm,4.0)-1.0);
				//printf("derivs: %d %lg %lg %lg %lg\n", i, G.CP[i], G.LoverT[i]*T[i], dxdT, G.LoverT[i]*T[i]*dxdT);
				G.CP[i]+= G.LoverT[i]*T[i] * dxdT;
			}
		}
	}
/*
	if (T[i] > G.Tmelt[i] && T[i+1] < G.Tmelt[i+1]) {
		double latentHeat = 0.805 * 1.38e-16*T[i]/(EOS.A[1]*1.67e-24);
	//	dTdt[i] += latentHeat/(ODE.timestep*G.CP[i]);
	dTdt[i] += -(latentHeat * 4.0 * G.y[i] *dTdt[i]/T[i])/(G.dx*G.CP[i]*G.y[i]);
	//	printf("here2 %d %lg %lg %lg %lg %lg %lg %lg\n", i, G.y[i], T[i], G.Tmelt[i], T[i-1], G.Tmelt[i-1],latentHeat,0.0);
}
*/
		
	for (int i=1; i<=G.N; i++) {
  		dTdt[i]=G.g*(G.F[i+1]-G.F[i])/(G.dx*G.CP[i]*G.P[i]);
		if (G.nuflag) dTdt[i]+=-(G.NU[i]/G.CP[i]);
		if (G.accreting) dTdt[i]+=G.EPS[i]/G.CP[i];
	}
  	dTdt[G.N+1]=0.0;

}

double calculate_heat_flux(int i, double *T)
{
	double flux;
	if (i>1 || (G.accreting && G.outburst_duration > 0.0))   
//		if (i>1 || (G.accreting && EOS.B == 0.0))   
		// use this inside the grid, or at the surface when we are accreting (which 
		// fixes the outer temperature)
	 	flux = 0.5*(G.K[i]+G.K[i-1])*(T[i]-T[i-1])/G.dx;	
	else {
		// cooling boundary condition
		if (EOS.B == 0.0) {
			// from my envelope calculation (makegrid.cc)
			flux = (G.g/2.28e14)*TEFF.get(T[i]);
		} else {
			// for magnetars we use
			// Potekhin & Yakovlev 2001 eq.(27)
			double T9 = T[i]*1e-9;
			double xi = T9 - 0.001*pow(1e-14*G.g,0.25)*sqrt(7.0*T9);
			flux = 5.67e-5 * 1e24 * G.g*1e-14 * (pow(7*xi,2.25)+pow(0.333*xi,1.25));
		
			// or use makegrid.cc calculation
			//flux = (G.g/2.28e14)*TEFF.get(T[i]);
			
			
			// now correct for B ... 
			if (G.angle_mu >= 0.0) {
				// use the enhancement along the field direction
				double chi1,chi2;
				double B12=EOS.B*1e-12;
				chi1 = 1.0 + 0.0492*pow(B12,0.292)/pow(T9,0.24);
			//	chi2 = sqrt(1.0 + 0.1076*B12*pow(0.03+T9,-0.559))/
			//		pow(1.0+0.819*B12/(0.03+T9),0.6463);
				double fcond = 4.0*G.angle_mu*G.angle_mu/(1.0+3.0*G.angle_mu*G.angle_mu);		
				flux *= fcond*pow(chi1,4.0);//+(1.0-fcond)*pow(chi2,4.0);

			} else {
				// or use eq. (31) or PY2001  which gives F(B)/F(0)
				double fac, a1,a2,a3,beta;
				beta = 0.074*sqrt(1e-12*EOS.B)*pow(T9,-0.45);
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

double dTdt(int i, double *T)
// calculates the time derivative for grid cell i (used when calculating the jacobian)
{
	int k=i-1; if (k<1) k=1;
	int k2=i+1; if (k2>G.N+1) k2=G.N+1;
	for (int j=k; j<=k2; j++) calculate_vars(j,T[j], G.P[j], &G.CP[j], &G.K[j], &G.NU[j],&G.EPS[j]);
	if (i==1) outer_boundary(T[1],G.K[1],G.CP[1],G.NU[1],G.EPS[1],&T[0],&G.K[0],&G.CP[0],&G.NU[0], &G.EPS[0]);
	if (i==G.N) inner_boundary(T[G.N],G.K[G.N],G.CP[G.N],G.NU[G.N], G.EPS[G.N],
										&T[G.N+1],&G.K[G.N+1],&G.CP[G.N+1],&G.NU[G.N+1],&G.EPS[G.N+1]);

	if (G.include_latent_heat) {
		double Gm = G.GammaT[i-1]/(0.5*(T[i-1]+T[i]));
		double Gp = G.GammaT[i]/(0.5*(T[i]+T[i+1]));
		if (Gm < 175.0 && Gp > 175.0) {
			double dxdT = (4.0/T[i]) * pow(175.0/Gm,4.0)/ (pow(Gp/Gm,4.0)-1.0);
			G.CP[i]+= G.LoverT[i]*T[i] * dxdT;
		}
	}

	double f=G.g*(calculate_heat_flux(i+1,T)-calculate_heat_flux(i,T))/(G.dx*G.CP[i]*G.P[i]);
	if (G.nuflag) f+=-(G.NU[i]/G.CP[i]);
	if (G.accreting) f+=G.EPS[i]/G.CP[i];
	
	/*
	if (G.include_latent_heat) {
		if (T[i] > G.Tmelt[i] && T[i+1] < G.Tmelt[i+1]) {
			double latentHeat = 0.805 * 1.38e-16*T[i]/(EOS.A[1]*1.67e-24);
			double df=-(latentHeat * 4.0 * G.y[i] *f/T[i]);
		//	printf("here %d %lg %lg %lg %lg %lg %lg %lg\n", i, G.y[i], T[i], G.Tmelt[i], T[i-1], G.Tmelt[i-1], latentHeat,df);
			f += df/(G.dx*G.CP[i]*G.y[i]);
			//f += latentHeat/(ODE.timestep*G.CP[i]);
			
		}
		}
		*/
	return f;
}

void outer_boundary(double T1, double K1, double CP1, double NU1, double EPS1,
		double *T0, double *K0, double *CP0, double *NU0, double *EPS0)  
{
	if (G.accreting && G.Tt>0.0) *T0=G.Tt;   // constant temperature during accretion
	else *T0=T1*(8.0-G.dx)/(8.0+G.dx);   // flux propto T^4
	*K0=K1; *CP0=CP1;
	if (G.nuflag) *NU0=NU1; else *NU0=0.0;
	if (G.accreting) *EPS0=EPS1; else *EPS0=0.0;
}

void inner_boundary(double TN, double KN, double CPN, double NUN, double EPSN,
	double *TN1, double *KN1, double *CPN1, double *NUN1, double *EPSN1)
{
	*TN1=G.Tc;   // fixed core temperature
//	*TN1=TN;    // zero flux inner boundary
	*KN1=KN;
	*CPN1=CPN;  
	if (G.nuflag) *NUN1=NUN; else *NUN1=0.0;
	if (G.accreting) *EPSN1=EPSN; else *NUN1=0.0;
	
}


void jacobn(double t, double *T, double *dfdt, double **dfdT, int n)
// calculates the Jacobian numerically
{
  double f, e=0.001;
  // this assumes the arrays dfdt and dfdT are preinitialized to zero (I changed odeint to do this)
  for (int i=1; i<n; i++) {
    if (i>1) {
      T[i-1]*=1.0+e; f=dTdt(i,T);
      T[i-1]/=1.0+e; dfdT[i][i-1]=(f-dfdt[i])/(T[i-1]*e);
    }
    T[i]*=1.0+e; f=dTdt(i,T);
    T[i]/=1.0+e; dfdT[i][i]=(f-dfdt[i])/(T[i]*e);
    if (i<=n) {
      T[i+1]*=1.0+e; f=dTdt(i,T);
      T[i+1]/=1.0+e; dfdT[i][i+1]=(f-dfdt[i])/(T[i+1]*e);
    }
  }
}  




void calculate_vars(int i, double T, double P, double *CP, double *K, double *NU,double *EPS)
{
	// sometimes we get a nan value for T here from the integrator
	// In this case, set the temperature to be some value.. this seems to
	// deal with this problem ok
	if (isnan(T)) T=1e7;
	
	double beta=log10(T);
	if (beta <= G.betamax && beta >= G.betamin) {
		// lookup values in the precalculated table
		int j = 1 + (int) ((beta-G.betamin)/G.deltabeta);
		double interpfac=(beta-(G.betamin + (j-1)*G.deltabeta))/G.deltabeta;
		// interpolate the thermal conductivity to the current
		// value of impurity parameter Q
		double K0,K1;//,K0perp,K1perp;
		K0=G.K0_grid[i][j] + (G.K0_grid[i][j+1]-G.K0_grid[i][j])*interpfac;
		K1=G.K1_grid[i][j] + (G.K1_grid[i][j+1]-G.K1_grid[i][j])*interpfac;
		//K0perp=G.K0perp_grid[i][j] + (G.K0perp_grid[i][j+1]-G.K0perp_grid[i][j])*interpfac;
		//K1perp=G.K1perp_grid[i][j] + (G.K1perp_grid[i][j+1]-G.K1perp_grid[i][j])*interpfac;
		// use something like this next line to hardwire Q values
		double Qval;
		if (G.hardwireQ) {
			if (G.rho[i] > G.Qrho) Qval=G.Qinner; else Qval=EOS.Q;
//			if (P>2.28e29) Qval=G.Qinner; else Qval=EOS.Q;
		} else {
			Qval = G.Qimp[i];	
		}
		double KK,KKperp;
		KK=G.g*K0*K1/(K0*Qval+(1.0-Qval)*K1);
		if (EOS.B > 0) {
		//KKperp = G.g*K0perp*K1perp/(K0perp*Qval+(1.0-Qval)*K1perp);		
		double fcond = 4.0*G.angle_mu*G.angle_mu/(1.0+3.0*G.angle_mu*G.angle_mu);		
		*K=fcond*KK;//+(1.0-fcond)*KKperp;	
	} else {
		*K=KK;
	}
		
		*CP=G.CP_grid[i][j] + (G.CP_grid[i][j+1]-G.CP_grid[i][j])*interpfac;
		if (G.nuflag) *NU=G.NU_grid[i][j] + (G.NU_grid[i][j+1]-G.NU_grid[i][j])*interpfac; 
		else *NU=0.0;
		if (G.accreting) {
			*EPS=G.EPS_grid[i][1];  // assume heating is independent of temperature 
		//	*EPS=(G.EPS_grid[i][j] + (G.EPS_grid[i][j+1]-G.EPS_grid[i][j])*interpfac); 
			*EPS=*EPS * G.mdot * G.g;
		}
		else *EPS=0.0;
	} else {
		// if beta is outside the range of the precalculated table, then calculate directly
		//printf("Note: log10T outside range (%lg)\n", beta);
		EOS.P=P; 
		EOS.T8=1e-8*T; 		
		EOS.rho=G.rho[i];
		set_composition();
		*CP=EOS.CP();
		
		double Qval=EOS.Q;
		if (G.hardwireQ) {
			if (G.rho[i] > G.Qrho) Qval=G.Qinner;
		} else {
			Qval = G.Qimp[i];	
		}
			
		double Qkeep = EOS.Q;
		EOS.Q=Qval;
		//if (EOS.B == 0.0) { 
		//	*K=EOS.rho*EOS.K_cond(0.0)*G.g/P;	
	//	} else {
			double Kcond, Kcondperp;
			Kcond = potek_cond(&Kcondperp);	
			double fcond = 4.0*G.angle_mu*G.angle_mu/(1.0+3.0*G.angle_mu*G.angle_mu);		
			Kcond = fcond*Kcond;// + (1.0-fcond)*Kcondperp;
			*K=EOS.rho*Kcond*G.g/P;
	//	}
		EOS.Q=Qkeep;
			
		//*K=3.03e20*pow(EOS.T8,3)/(EOS.opac()*y);
		if (G.nuflag) *NU=EOS.eps_nu(); else *NU=0.0;	
		if (G.accreting) *EPS=crust_heating(i)*energy_deposited(i)*G.mdot * G.g; else *EPS=0.0;
	}
 }


// --------------------------------------------------------------------------



// --------------------------------------------------------------------------
// Routines for initial setup 

void set_up_initial_temperature_profile_piecewise(char *fname)
{
	// first read in the specified temperature-density relation from the file
	FILE *fp;
	char s1[100];
	printf("Reading initial temperature profile from %s\n",fname);
	fp = fopen(fname,"r");
	double *rhovec, *Tvec;
	rhovec=vector(1,100);
	Tvec=vector(1,100);
	int i=1;
	while (!feof(fp)) {		
		double rho, T;
		// new lines: read from lines marked ">" in init.dat
		char *e=fgets(s1,200,fp);		
		if (!strncmp(s1,">",1)) {
			sscanf(s1,">%lg\t%lg\n",&rho,&T);
			// old: direct read from Tinit.dat
			//		fscanf(fp,"%lg %lg\n",&rho,&T);
			if (rho>=-1.0) {
				if (rho == 0.0) rho=G.rho[1];
				if (rho == -1.0) rho=G.rho[G.N];
				rhovec[i] = rho;
				if (T < 0.0) T=G.Tc;
				Tvec[i] = T;
				i++;
			}
		}
	}
	fclose(fp);
	int nvec=i-1;

	double totalEd=0.0;

	// now assign initial temperatures
	if (G.output) fp = fopen("gon_out/initial_condition","w");
	double I=0.0;   // I is the integral used to get the thermal time
	double E=0.0;   // total energy deposited
	for (int i=1; i<=G.N+1; i++) {

		double Ti;
		if (i==1) {
			Ti=Tvec[1]; 
			G.Tt=Ti;
		} else {
			if (i==G.N+1) {
				Ti=Tvec[nvec];		
				//G.Tc=Ti;
			} else {
				int	j=1; 
				while (rhovec[j] < G.rho[i]) j++;
				Ti = pow(10.0,log10(Tvec[j-1]) + log10(Tvec[j]/Tvec[j-1])*log10(G.rho[i]/rhovec[j-1])/log10(rhovec[j]/rhovec[j-1]));
				//printf("%lg %lg %lg %lg %lg\n",Tvec[j-1],Tvec[j],rhovec[j-1],rhovec[j],G.rho[i]);
			}
		}	
		
		ODE.set_bc(i,Ti);							
				
		// Figure out the energy deposited
		EOS.P=G.P[i];
		EOS.T8=G.Tc*1e-8;
		EOS.rho=EOS.find_rho();
		set_composition();				
		Ode_Int ODEheat;
		ODEheat.init(1);
		ODEheat.set_bc(1,0.0);
		double Ed;
		ODEheat.go(G.Tc,Ti,0.01*G.Tc, 1e-6, heatderivs2);
		Ed = ODEheat.get_y(1,ODEheat.kount);
		totalEd+=Ed*4.0*PI*1e10*G.radius*G.radius*G.dx*G.P[i]/G.g;
		if (Ed > 0.0) printf("heating cell %d:  Ti=%lg Tf=%lg E25=%lg rho=%lg\n",
					i, G.Tc, Ti, Ed*G.rho[i]*1e-25, G.rho[i]);
		ODEheat.tidy();			
					
		if (G.output) {
			EOS.P=G.P[i]; EOS.rho=G.rho[i]; EOS.T8=1e-8*Ti;
			set_composition();				

			double TTC=EOS.TC();
			if (EOS.Yn >0.0 && Ti > TTC && TTC > Tvec[nvec]) {
				printf("%d %lg %lg %lg\n", i, TTC, Ti, EOS.Yn);
				ODE.set_bc(i,TTC);
			}

			double null=0.0;
			double Kcond = potek_cond(&null);

			I+=sqrt(EOS.CP()/(Kcond*EOS.rho))*(G.P[i]-G.P[i-1])/G.g;
			double tt = I*I*0.25/(24.0*3600.0);
			
			fprintf(fp, "%d %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg\n", i, G.P[i], Ti, EOS.rho,EOS.CV(), 
					Kcond, EOS.Yn, 1e-39*EOS.Yn * EOS.rho/1.67e-24, tt,E,EOS.A[1],EOS.Z[1],TTC);
		}
	}	
	
	printf("Total energy input to get this initial T profile = %lg\n",totalEd);
	
	if (G.output) fclose(fp);
	
	free_vector(rhovec,1,100);
	free_vector(Tvec,1,100);
	
}





void set_up_initial_temperature_profile(void)
// initialize the temperature profile on the grid
{
	
	// first run some accretion with heating, long enough for the crust
	// to get into a thermal steady-state
	G.accreting = !G.instant_heat;
	for (int i=G.N+1; i>=1; i--) {
		// a linear profile between top and bottom
		//double Ti = pow(10.0,log10(G.Tc) + log10(G.Tt/G.Tc)*log10(G.P[i]/G.Pb)/log10(G.Pt/G.Pb));
		// or constant profile
		double Ti = G.Tc;
		// a linear profile adjusts to steady state *much* more quickly,
		// but for XTEJ for example I want to heat up from isothermal and the crust 
		// does not get to steady state
		ODE.set_bc(i,Ti);
	}
	clock_t timer;
	start_timing(&timer);
	double dt=G.outburst_duration*3.15e7*0.01;
	if (dt > 1e6) dt=1e6;
	ODE.go(0.0, G.outburst_duration*3.15e7,dt, 1e-8, derivs);
	stop_timing(&timer,"ODE.go (initial heating)");
	printf("number of steps = %d  (time=%lg)\n", ODE.kount, ODE.get_x(ODE.kount));

	// set initial condition and write out some info
	FILE *fp;
	if (G.output) fp = fopen("gon_out/initial_condition","w");
	double I=0.0;   // I is the integral used to get the thermal time
	double E=0.0;   // total energy deposited
	for (int i=1; i<=G.N+1; i++) {
		
		double Ti = ODE.get_y(i,ODE.kount);
		if (crust_heating(i) > 0.0 && G.instant_heat) {				
			// find the initial temperature from instantaneous heating
			EOS.P=G.P[i];
			set_composition();				
			Ode_Int ODEheat;
			ODEheat.init(1);
			ODEheat.set_bc(1,Ti);
			printf("heating cell %d:  Ti=%lg E25=%lg ", i,Ti,energy_deposited(i));
			ODEheat.go(0.0, energy_deposited(i), 1e-4, 1e-6, heatderivs);
			Ti = ODEheat.get_y(1,ODEheat.kount);
			printf(" Tf=%lg\n", Ti);
			ODEheat.tidy();			
		}
		
		ODE.set_bc(i,Ti);
					
		if (G.output) {
			EOS.P=G.P[i]; EOS.rho=G.rho[i]; EOS.T8=1e-8*Ti;
			set_composition();				
				//double f=1.0+EOS.Uex()/3.0;
				//double Pion=8.254e15*this->rho*this->T8*Yi()*f;   // ions
			   //double Pe=EOS.pe();
					
			double Qval;			
				if (G.hardwireQ) {
					if (G.rho[i] > G.Qrho) Qval=G.Qinner; else Qval=EOS.Q;
						//			if (P>2.28e29) Qval=G.Qinner; else Qval=EOS.Q;
								} else {
									Qval = G.Qimp[i];	
								}
								double Q_store=EOS.Q;
								EOS.Q=Qval;
								
								double Kcondperp=0.0;
			double Kcond = potek_cond(&Kcondperp);
								
			I+=sqrt(EOS.CP()/(Kcond*EOS.rho))*(G.P[i]-G.P[i-1])/G.g;
			double tt = I*I*0.25/(24.0*3600.0);
			
			E+=energy_deposited(i)*crust_heating(i)*G.mdot*G.g*G.outburst_duration*3.15e7*
				4.0*PI*G.radius*G.radius*1e10*(G.P[i]-G.P[i-1])/G.g;
				
			fprintf(fp, "%d %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg\n", i, G.P[i], Ti, EOS.rho,EOS.CV(), 
					Kcond, EOS.Yn, 1e-39*EOS.Yn * EOS.rho/1.67e-24, tt,E,EOS.A[1],EOS.Z[1],EOS.TC(),
					EOS.Chabrier_EF(),Kcondperp);
				EOS.Q=Q_store;	
		}
	}	
	if (G.output) fclose(fp);

	printf("Total energy deposited=%lg\n",E);

}







void precalculate_vars(void) 
// calculate various quantities at each grid point as a function of temperature
// then during the run we can look them up in a table
{
	// the table is constructed in terms of log10(T)
	// for historical reasons, this is called beta here
	// (for long X-ray bursts where radiation pressure is significant,
	// beta=Prad/P is a better variable to use)
	G.nbeta=100;
	G.betamin=7.0;
	G.betamax=10.0;
	G.deltabeta = (G.betamax-G.betamin)/(1.0*(G.nbeta-1));	

	G.CP_grid = matrix(1,G.N+2,1,G.nbeta);	
	G.K1_grid = matrix(1,G.N+2,1,G.nbeta);	
	G.K0_grid = matrix(1,G.N+2,1,G.nbeta);	
	G.K1perp_grid = matrix(1,G.N+2,1,G.nbeta);	
	G.K0perp_grid = matrix(1,G.N+2,1,G.nbeta);	
	G.NU_grid = matrix(1,G.N+2,1,G.nbeta);	
	G.EPS_grid = matrix(1,G.N+2,1,G.nbeta);	

	char s[100];
	if (EOS.B > 0.0) sprintf(s,"gon_out/precalc_results_%lg",log10(EOS.B));
	else sprintf(s,"gon_out/precalc_results_0");

	int read_grid=0;
	FILE *fp;
	fp=fopen(s,"r");
	if (fp == NULL || G.force_precalc ==1) {
		fp=fopen(s,"w");
		printf("Precalculating quantities and writing to file %s...\n",s);

		for (int i=1; i<=G.N+1; i++) {
	
			EOS.P=G.P[i];
			EOS.rho = G.rho[i];
			set_composition();
			
			fprintf(fp, "Grid point %d  P=%lg  rho=%lg  A=%lg  Z=%lg Yn=%lg:  T8,CP,K,eps_nu,eps_nuc\n",
				i, G.P[i], G.rho[i], (1.0-EOS.Yn)*EOS.A[1], EOS.Z[1], EOS.Yn);
		
			double heating = crust_heating(i);
		
			for (int j=1; j<=G.nbeta; j++) {		
				double beta = G.betamin + (j-1)*(G.betamax-G.betamin)/(1.0*(G.nbeta-1));
				EOS.T8 = 1e-8*pow(10.0,beta);

				G.CP_grid[i][j]=EOS.CV();
				G.NU_grid[i][j]=EOS.eps_nu();			
				G.EPS_grid[i][j]=heating;

				// we calculate the thermal conductivity for Q=0 and Q=1, and later interpolate to the
				// current value of Q. This means we can keep the performance of table lookup even when
				// doing MCMC trials which vary Q.
				double Q_store=EOS.Q;  // store Q temporarily
				EOS.Q=0.0;
				double Kcond,Kcondperp;
				//Kcond = EOS.K_cond(EOS.Chabrier_EF());
				Kcond = potek_cond(&Kcondperp);   
				Kcondperp=0.0;
				G.K0_grid[i][j]=EOS.rho*Kcond/G.P[i];
				G.K0perp_grid[i][j]=EOS.rho*Kcondperp/G.P[i];
				EOS.Q=1.0;
				//Kcond = EOS.K_cond(EOS.Chabrier_EF());
				Kcond = potek_cond(&Kcondperp);
				Kcondperp=0.0;
				G.K1_grid[i][j]=EOS.rho*Kcond/G.P[i];
				G.K1perp_grid[i][j]=EOS.rho*Kcondperp/G.P[i];
				EOS.Q=Q_store;  // restore to previous value

				// We include thermal conductivity only. The equivalent expression but putting in all sources
				// of opacity would be
				//3.03e20*pow(EOS.T8,3)/(EOS.opac()*G.y[i]);

				fprintf(fp, "%lg %lg %lg %lg %lg %lg %lg %lg\n", EOS.T8, G.CP_grid[i][j], 
					G.K0_grid[i][j],G.K1_grid[i][j], G.K0perp_grid[i][j],G.K1perp_grid[i][j],
					G.NU_grid[i][j], G.EPS_grid[i][j]);

				G.EPS_grid[i][j]*=energy_deposited(i);

			}	
		}

	} else {
		
		printf("***Reading precalculated quantities from file %s\n", s);
		
		for (int i=1; i<=G.N+1; i++) {
			int kk; double dd;
			fscanf(fp, "Grid point %d  P=%lg  rho=%lg  A=%lg  Z=%lg Yn=%lg:  T8,CP,K,eps_nu,eps_nuc\n",
					&kk,&dd,&dd,&dd,&dd,&dd);
			for (int j=1; j<=G.nbeta; j++) {		
				fscanf(fp, "%lg %lg %lg %lg %lg %lg %lg %lg\n", &EOS.T8, &G.CP_grid[i][j], 
					&G.K0_grid[i][j],&G.K1_grid[i][j], &G.K0perp_grid[i][j],&G.K1perp_grid[i][j],
					&G.NU_grid[i][j], &G.EPS_grid[i][j]);
				// always calculate the crust heating..
				G.EPS_grid[i][j]=crust_heating(i);
				G.EPS_grid[i][j]*=energy_deposited(i);
			}
		}
			
	}
	
	fclose(fp);
}


double energy_deposited(int i)
{
	if (G.rho[i]>4e11) return G.energy_deposited_inner;
	else return G.energy_deposited_outer;
}


double crust_heating(int i) 
// calculates the crust heating in erg/g/s
// for grid point i
{
	double eps=0.0;

	// if B>0 we are doing a magnetar
	// changed this to:
	// if we are heating on < 1day timescale then its a magnetar
	if (G.outburst_duration<1.0/365.0) {
		
		// eps in erg/g/s
		eps = 1e25/(G.rho[i]*G.outburst_duration*3.15e7);
		eps /= G.mdot * G.g;   // modify to the units used in the code
		
//		if (G.rho[i]>1e9) eps *= pow(G.rho[i]/1e9,0.4);
//		eps *= G.rho[i]/1e11;
		
		// the next line limits the heating to a limited region of the crust
		if (G.rho[i] < G.rhot) eps=0.0;
		if (G.rho[i] > G.rhob) eps=0.0;
		
	} else {  // otherwise we are doing an accreting neutron star

		if (!G.hardwireQ) 
			eps = G.Qheat[i]*8.8e4*9.64e17/(G.P[i]*G.dx);
		else
			eps = crust_heating_rate(G.P[i]);
	}

	return eps;	
}


double crust_heating_rate(double P) 
// calculates the crust heating in erg/g/s
 // except for a factor of gravity --- multiply by gravity to get these units
{
	double eps=0.0;
	// simple "smeared out" heating function, 1.2MeV in inner crust, 0.2MeV in outer crust
	if (P >= 6e15*2.28e14 && P <= 1e17*2.28e14) eps=8.8e4*1.2*9.64e17/(P*log(1e17/6e15));
// 	if (y >= 1e12 && y < 1e15) eps=G.mdot*8.8e4*0.15*9.64e17/(y*log(1e15/1e12));
 	if (P >= 3e12*2.28e14 && P < 3e15*2.28e14) eps=8.8e4*0.2*9.64e17/(P*log(3e15/3e12));
//if (y >= 6e15 && y <= 3e18) eps=G.mdot*8.8e4*1.2*9.64e17/(y*log(3e18/6e15));

	return eps;	
}



void get_TbTeff_relation(void)
// reads in the Flux-T relation from the data file output by makegrid.cc
{
	double *temp, *flux;  // temporary storage to initialize the spline
	int npoints = 195;  //  needs to be >= number of points read in
	temp = vector(1,npoints);
	flux = vector(1,npoints);
	
	// the file "out/grid" is made by makegrid.cc
	// it contains  (column depth, T, flux)  in cgs
	FILE *fp;
	if (G.gpe) fp = fopen("out/grid_He4","r");
	else fp = fopen("out/grid_He9","r");
	//FILE *fp = fopen("grid_sorty","r");
	FILE *fp2 = fopen("gon_out/TbTeff", "w");
	
	double y,T,F,rho;
	int count = 0;
	while (!feof(fp)) {
		
		fscanf(fp, "%lg %lg %lg %lg\n", &y, &T, &F,&rho);
		//printf("%lg %lg %lg\n", y, T, F);
		if (y == log10(G.yt)) {  // select out the points which correspond to the top column
			count++;
			temp[count] = pow(10.0,T);
			// correct for gravity here:
			flux[count] = pow(10.0,F); //* (G.g/2.28e14);
			//printf("%lg %lg\n", T, F);
			fprintf(fp2, "%d %lg %lg %lg %lg %lg\n", count,y,T,F,temp[count],flux[count]);
		}
	}
		
	fclose(fp);
	fclose(fp2);
	
	// the following spline contains the flux as a function of temperature
	// at column depth G.yt
	TEFF.minit(temp,flux,count);
	
	free_vector(temp,1,npoints);
	free_vector(flux,1,npoints);
}



void set_up_grid(int ngrid, const char *fname)
// allocates storage for the grid and also computes the density
// and composition at each grid point
{
  	G.N=ngrid;   // number of grid points
	G.Pb=6.5e32; // column depth at the crust/core boundary
						  // we used to set this to y=3e18 but now fix pressure
  	G.Pt=G.yt*2.28e14;   // pressure at the top

	Spline QiSpline;
	Spline QhSpline;
	{
		FILE *fp=fopen(fname,"r");
	
		int npoints;
		fscanf(fp,"%d",&npoints);
		npoints--;
		printf("Crust model has %d points\n",npoints);
		
		double *Qi,*Qh,*AA,*ZZ,*rho,*Yn;
		Qi=vector(1,npoints);
		Qh=vector(1,npoints);
		AA=vector(1,npoints);
		ZZ=vector(1,npoints);
		Yn=vector(1,npoints);
		rho=vector(1,npoints);

		for (int i=1; i<=npoints; i++) {
			double dummy;
			fscanf(fp, "%lg %lg %lg %lg %lg %lg %lg %lg\n",
				&dummy,&rho[i],&dummy,&Qh[i],&ZZ[i],&AA[i],&Qi[i],&Yn[i]);
			//printf("%lg %lg %lg %lg %lg\n", rho[i], Qh[i],ZZ[i],AA[i],Qi[i]);
			rho[i]=log10(rho[i]);
		}

		YnSpline.minit(rho,Yn,npoints);
		ZZSpline.minit(rho,ZZ,npoints);
		AASpline.minit(rho,AA,npoints);
		QiSpline.minit(rho,Qi,npoints);
		QhSpline.minit(rho,Qh,npoints);
	/*	
		for (int i=0; i<=G.N+1; i++) {
			double x=log(G.yt)+G.dx*(i-1);
			double yh = exp(x+0.5*G.dx);
			double rhoh = pow(10.0,RHO.get(log10(G.g*yh)));
			EOS.rho = rhoh;
			set_composition();
			// GammaT[i] refers to i+1/2
			G.GammaT[i] = pow(EOS.Z[1]*4.8023e-10,2.0)*pow(4.0*PI*rhoh/(3.0*EOS.A[1]*1.67e-24),1.0/3.0)/1.38e-16;
		}
	*/	
		free_vector(Qi,1,npoints);
		free_vector(Qh,1,npoints);
		free_vector(AA,1,npoints);
		free_vector(ZZ,1,npoints);		
		free_vector(Yn,1,npoints);		
		fclose(fp);

	}



  	// We will compute density at each grid point
  	// the only complication is that composition depends on the local
  	// density, so the approach we take is to 
  	// first set up a spline P(rho) and then use this to 
  	// interpolate the correct density for each pressure
  	EOS.T8=1.0;  // have to set the temperature to be something
  	double *dens = vector(1,1001);
  	double *pres = vector(1,1001);
  	for (int i=1; i<=1001; i++) {
		dens[i] = 7.0 + (i-1)*(14.3-7.0)*0.001;
		EOS.rho = pow(10.0,dens[i]);
		set_composition();
		pres[i] = log10(EOS.ptot());
  	}
  	RHO.minit(pres,dens,1001);
  	free_vector(dens,1,1001);
  	free_vector(pres,1,1001);	

  	// storage
  	G.rho=vector(0,G.N+2);  
  	G.Tmelt=vector(0,G.N+2);  
  	G.CP=vector(0,G.N+2);
  	G.P=vector(0,G.N+2);
  	G.K=vector(0,G.N+2);
  	G.F=vector(0,G.N+2);
  	G.NU=vector(0,G.N+2);
  	G.EPS=vector(0,G.N+2);
  	G.GammaT=vector(0,G.N+2);
  	G.LoverT=vector(0,G.N+2);

  	G.Qheat=vector(0,G.N+2);
  	G.Qimp=vector(0,G.N+2);

  	G.dx=log(G.Pb/G.Pt)/(G.N-1);   // the grid is equally spaced in log column
  

	FILE *fp = fopen("gon_out/grid_profile","w");

	// calculate Gamma x T at the half grid points
	for (int i=0; i<=G.N+1; i++) {
		double x=log(G.Pt)+G.dx*(i-1);
		double Ph = exp(x+0.5*G.dx);
		double rhoh = pow(10.0,RHO.get(log10(Ph)));
		EOS.rho = rhoh;
		set_composition();
		// GammaT[i] refers to i+1/2
		G.GammaT[i] = pow(EOS.Z[1]*4.8023e-10,2.0)*pow(4.0*PI*rhoh/(3.0*EOS.A[1]*1.67e-24),1.0/3.0)/1.38e-16;
	}

	double Qtot=0.0;
  	for (int i=0; i<=G.N+2; i++) {
    	double x=log(G.Pt)+G.dx*(i-1);
    	G.P[i]=exp(x);

		double P1 = exp(x-0.5*G.dx);
		double P2 = exp(x+0.5*G.dx);
		double rho1 = pow(10.0,RHO.get(log10(P1)));
		double rho2 = pow(10.0,RHO.get(log10(P2)));
		G.Qheat[i]=QhSpline.get(log10(rho2))-QhSpline.get(log10(rho1));
		if (G.Qheat[i]<0.0) G.Qheat[i]=0.0;
		Qtot+=G.Qheat[i];
	
    	G.rho[i] = pow(10.0,RHO.get(log10(G.P[i])));
		EOS.rho = G.rho[i];
		set_composition();
		
		G.Tmelt[i] = 5e8*pow(G.P[i]/(2.28e14*1.9e13),0.25)*pow(EOS.Z[1]/30.0,5.0/3.0);
		G.LoverT[i] = 0.8 * 1.38e-16 /(EOS.A[1]*1.67e-24);
		if (EOS.Q>=0.0) {
			G.hardwireQ=1;
			// the Q values are assigned directly in 'calculate_vars'
			if (i==0) printf("Using supplied Qimp values and HZ composition and heating.\n");
		} else {
			G.hardwireQ=0;
			G.Qimp[i]=QiSpline.get(log10(G.rho[i]));
			if (G.Qimp[i] < 0.0) G.Qimp[i]=0.0;
			if (i==0) printf("Using Qimp, composition, and heating from the crust model.\n");
		}
		fprintf(fp, "%d %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg\n", i, G.P[i], G.rho[i], EOS.A[1]*(1.0-EOS.Yn), 
						EOS.Z[1], EOS.Yn,EOS.A[1],EOS.ptot(), G.Tmelt[i], G.GammaT[i]/1e8, G.LoverT[i]*1e8,
						AASpline.get(log10(G.rho[i])),  ZZSpline.get(log10(G.rho[i])), G.Qimp[i], G.Qheat[i]);
  	}

	fclose(fp);

	QiSpline.tidy();

  	printf("Grid has %d points, delx=%lg, Pb=%lg, rhob=%lg, Pt=%lg, rhot=%lg\n", 
			G.N, G.dx, G.P[G.N],G.rho[G.N],G.P[1],G.rho[1]);
		printf("Total heat release is %lg MeV\n",Qtot);
}


void set_composition(void)
{
	if (G.hardwireQ) {
		// use the EOS routines to get the composition
		// ie. crust models from the literature
		EOS.set_comp();	
	} else {
		// otherwise use our crust model
		// the model gives the mean A, mean Z and Yn
		// and we set up the variables as in our EOS.set_comp() routine
		EOS.Yn = YnSpline.get(log10(EOS.rho));
		if (EOS.Yn < 1e-6) EOS.Yn=0.0;
		EOS.A[1]=AASpline.get(log10(EOS.rho));
		EOS.A[1]/=(1.0-EOS.Yn);
		EOS.Z[1]=ZZSpline.get(log10(EOS.rho));
		EOS.set_Ye=EOS.Z[1]/EOS.A[1];		
		//printf("%lg %lg %lg %lg\n", EOS.Yn, EOS.rho, EOS.A[1], EOS.Z[1]);
	}
}

void set_ns_parameters(double mass, double radius)
	// mass in solar masses; radius in km
{
	radius*=1e5;
	G.ZZ=1.0/sqrt(1.0-2.0*6.67e-8*2e33*mass/(9e20*radius));
	G.g=G.ZZ*6.67e-8*2e33*mass/(radius*radius);	
	printf("NS parameters: M %lg M_sun, R %lg km, g14=%lg 1+z=%lg\n",
		mass, radius/1e5, G.g/1e14, G.ZZ);
}

void set_ns_redshift(double R)
// given gravity (G.g)  and radius in km
// sets the NS redshift
{
	double y, x;
	R*=1e5;
	y = pow(R*G.g/9e20,2.0);
	// x is  GM/Rc^2
	x = y*(sqrt(1.0+y)-1.0);
	G.mass = R*9e20*x/6.67e-8;
	G.mass/=2e33;
	G.ZZ=1.0/sqrt(1.0-2.0*x);	
}



//---------------------------------------------------------------------------------------
// wrapper for Potekhin's conductivity routine
double potek_cond(double *Kperp)
  {
      double s1,s2,s3,k1,k2,k3;
      double Bfield=EOS.B/4.414e13, Zimp=sqrt(EOS.Q), AA=EOS.A[1]*(1.0-EOS.Yn);
      double temp=EOS.T8*1e2/5930.0;
      double rr=EOS.rho/(EOS.A[1]*15819.4*1822.9);
      condegin_(&temp,&rr,&Bfield,&EOS.Z[1],&AA,&EOS.A[1],&Zimp, &s1,&s2,&s3,&k1,&k2,&k3);
      //return k2*2.778e15;
	double kk=k1*2.778e15;
	*Kperp = k2*2.778e15;
	//  k1, k2 and k3 are longitudinal, transverse and off-diagonal conductivities

//	return kk;

	double ksph =0.0;
	if (EOS.Yn > 0.0 && G.include_sph) {
		// now SF conductivity
		double lsph = 0.001 * pow(EOS.rho/4e11,-1.0/3.0);  // mfp in cm
		
		double vs = 1.05e-27/(sqrt(3.0)*1.67e-24*3e10);
		vs*=pow(3.0*PI*PI*EOS.rho*EOS.Yn/1.67e-24,1.0/3.0);	
		//printf("%lg %lg %lg %lg\n", EOS.T8, EOS.rho, EOS.Yn, vs);
					
//		double vs = 0.1 * pow(EOS.rho/4e11,1.0/3.0);
		ksph = 1.5e22 * pow(EOS.T8,3.0) * pow(0.1/vs,2.0) * lsph;
	}
	
	//ksph=0.0;
	
	*Kperp += ksph;
	return kk + ksph;

  }



void output_some_EOS_info()
{

	if (0)
	{
		FILE *fp=fopen("neut2.dat","w");

		EOS.A[1]=56.0;
		EOS.Z[1]=26.0;

		EOS.B=1e13;
		EOS.rho=1e11;

		for (int i=0; i<=100; i++) {

			EOS.T8=pow(10.0,8.0 + 2.0*0.01*i)*1e-8;

			fprintf(fp, "%lg %lg", EOS.T8, EOS.rho*EOS.eps_nu());
			fprintf(fp, " %lg", EOS.rho*EOS.Q1);
			fprintf(fp, " %lg", EOS.rho*EOS.Q2);
			fprintf(fp, " %lg", EOS.rho*EOS.Q3);
			fprintf(fp, " %lg", EOS.rho*EOS.Q4);
			fprintf(fp, " %lg", EOS.rho*EOS.Q5);
			fprintf(fp, " %lg\n", 1e-10);

		}

		fclose(fp);

		}


		if (0)
		{
			FILE *fp=fopen("neut.dat","w");

			EOS.A[1]=56.0;
			EOS.Z[1]=26.0;

			EOS.B=1e13;

			for (int i=0; i<=100; i++) {

				EOS.T8=pow(10.0,8.0 + 2.0*0.01*i)*1e-8;

				EOS.rho=1e9;
				fprintf(fp, "%lg %lg", EOS.T8, EOS.rho*EOS.eps_nu());

				EOS.rho=3e9;
				fprintf(fp, " %lg", EOS.rho*EOS.eps_nu());

				EOS.rho=1e10;
				fprintf(fp, " %lg", EOS.rho*EOS.eps_nu());

				EOS.rho=3e10;
				fprintf(fp, " %lg", EOS.rho*EOS.eps_nu());

				EOS.rho=1e11;
				fprintf(fp, " %lg", EOS.rho*EOS.eps_nu());

				EOS.rho=3e11;
				fprintf(fp, " %lg", EOS.rho*EOS.eps_nu());

				EOS.rho=1e12;
				fprintf(fp, " %lg\n", EOS.rho*EOS.eps_nu());

			}

			fclose(fp);

			}


if (0)
{
	FILE *fp=fopen("CV.dat","w");

	EOS.A[1]=56.0;
	EOS.Z[1]=26.0;

	EOS.B=1e8;

	EOS.accr=0;

	for (int i=0; i<=100; i++) {

		EOS.T8=pow(10.0,8.0 + 2.0*0.01*i)*1e-8;

		EOS.rho=1e9; EOS.set_comp();
		
		fprintf(fp, "%lg %lg", EOS.T8, EOS.CP());

		EOS.rho=3e9; EOS.set_comp();
		fprintf(fp, " %lg", EOS.CP());

		EOS.rho=1e10; EOS.set_comp(); 
		fprintf(fp, " %lg", EOS.CP());

		EOS.rho=3e10; EOS.set_comp();
		fprintf(fp, " %lg", EOS.CP());

		EOS.rho=1e11; EOS.set_comp();
		fprintf(fp, " %lg", EOS.CP());

		EOS.rho=3e11; EOS.set_comp();
		fprintf(fp, " %lg", EOS.CP());

		EOS.rho=1e12; EOS.set_comp();
		fprintf(fp, " %lg\n",EOS.CP());

	}

	fclose(fp);

	}

	if (0)
	{
		FILE *fp=fopen("TC.dat","w");

	 	EOS.accr = 0;   // set crust composition
		EOS.T8=1.0;

		for (int i=0; i<=100; i++) {
			EOS.rho=pow(10.0,11.0 + 3.0*0.01*i);
			EOS.set_comp();

			EOS.gap=1;
			fprintf(fp, "%lg %lg", EOS.rho, EOS.TC());

			EOS.gap=5;
			fprintf(fp, " %lg", EOS.TC());

			EOS.gap=6;
			fprintf(fp, " %lg\n", EOS.TC());
		}

		fclose(fp);

	}
	
	
	
	if (0)
	{
		FILE *fp=fopen("heat.dat","w");

		EOS.B=1e15;
		EOS.accr=0;
		EOS.rho = 1e9;
		EOS.set_comp();

		EOS.T8=0.7;
		EOS.P = EOS.ptot();

		ODE.init(1);
		ODE.set_bc(1,7e7);
		ODE.go(0.0, 3000.0, 1e-4, 1e-6, heatderivs);
		
		for (int i=1; i<=ODE.kount; i++) {

			EOS.T8 = ODE.get_y(1,i)*1e-8;
			EOS.P = G.Pb; 
			EOS.rho = EOS.find_rho();
			fprintf(fp,"%lg %lg %lg\n",ODE.get_x(i),ODE.get_y(1,i),EOS.rho);

		}

		fclose(fp);

		printf("Finished heat integration.\n");

	}

	
	
}




void heatderivs2(double T, double E[], double dEdT[])
// calculates the time derivatives for the whole grid
{
	EOS.T8 = T/1e8;
	EOS.rho = EOS.find_rho();
	
//	printf("%lg %lg %lg %lg %lg %lg %lg\n", T, EOS.P, EOS.rho,E[1], EOS.A[1], EOS.Z[1], EOS.X[1]);
	
	dEdT[1] = EOS.CP();
}

void heatderivs(double E, double T[], double dTdE[])
// calculates the time derivatives for the whole grid
{
	EOS.T8 = T[1]/1e8;
	EOS.rho = EOS.find_rho();
	
	dTdE[1] = 1e25/(EOS.rho*EOS.CP());
}



void make_TbTeff_relation(void)
{
	Ode_Int SURFODE;
	FILE *fp;
	
	fp=fopen("out/mygrid","w");

	EOS.X[1]=1.0; EOS.A[1]=56.0; EOS.Z[1]=26.0;   // iron composition
	EOS.Yn=0.0; 
	//EOS.Q=0.0; EOS.B=0.0;

	G.surfy=-1.0;

	for (int j=0; j<=16; j++) {

		G.surfF = pow(10.0,17.0+0.5*j);
			
		SURFODE.init(1);
		
		//double yt;
	  	//yt=zbrent(find_surf_eqn,1e-6,10.0,1e-5);
	  	//if (yt==1e-6 || yt==10.0) printf("yt out of bounds (%lg)\n", yt);
		//printf("%lg %lg\n", G.surfF, yt);
	 		
		double yt = 0.0;
		double Touter = pow(G.surfF/5.67e-5,0.25);
		SURFODE.set_bc(1,Touter);	
		SURFODE.go_simple(yt,14.0,(int)((14.0-yt)/0.01),surfderivs);
		//SURFODE.go(0.0,14.0,0.1,1e-6,surfderivs);
		for (int i=1; i<=SURFODE.kount; i++) {

			EOS.T8 = SURFODE.get_y(1,i)*1e-8;
			EOS.P = G.g * pow(10.0,SURFODE.get_x(i)); 
			if (SURFODE.get_x(i) < G.surfy){
				EOS.A[1]=4.0; EOS.Z[1]=2.0;
			} else {
				EOS.A[1]=56.0; EOS.Z[1]=26.0;
			}
			EOS.rho = EOS.find_rho();
			double kapp = EOS.opac();
			double null = 0.0;
			kapp = 1.0/(EOS.kff + EOS.kes) + 1.0/(3.03e20*pow(EOS.T8,3.0)/(potek_cond(&null)*EOS.rho));
			kapp = 1.0/kapp;
			fprintf(fp,"%lg %lg %lg %lg %lg %lg\n",
				SURFODE.get_x(i),log10(SURFODE.get_y(1,i)),log10(G.surfF),log10(EOS.rho),
				EOS.opac(), kapp);

		}
		SURFODE.tidy();
	}

	fclose(fp);
	printf("Made Tb-Teff grid.\n");

}



double find_surf_eqn(double y)
{
	EOS.P=G.g*y; EOS.T8=1e-8*pow(G.surfF/5.67e-5,0.25);
	EOS.rho=EOS.find_rho();
	return EOS.opac()*y-2.0/3.0;
}


void surfderivs(double y, double T[], double dTdy[])
// calculates the time derivatives for the whole grid
{
	EOS.T8 = T[1]/1e8;
	EOS.P=G.g*pow(10.0,y);
	if (y < G.surfy) {
		EOS.A[1]=4.0; EOS.Z[1]=2.0;
	} else {
		EOS.A[1]=56.0; EOS.Z[1]=26.0;
	}
	EOS.rho = EOS.find_rho();

	double kapp = EOS.opac();
	//kapp = 1.0/(EOS.kff + EOS.kes) + 1.0/(3.03e20*pow(EOS.T8,3.0)/(potek_cond()*EOS.rho));
	//kapp = 1.0/kapp;
	
	double norm = 16.0*5.67e-5*1e24*pow(EOS.T8,3.0)/(3.0*kapp);	
	dTdy[1] = 2.303*pow(10.0,y)*G.surfF/norm;
}


