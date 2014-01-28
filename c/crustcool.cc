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
double crust_heating_rate(int i);
double crust_heating(int i);
double energy_deposited(int i);
void output_result_for_step(int j, FILE *fp, FILE *fp2, double timesofar, double *last_time_output);
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
void heatderivs(double t, double T[], double dTdt[]);
void heatderivs2(double t, double T[], double dTdt[]);

// global variables
struct globals {
  	int N;  
  	double dx; 
  	double *P, *CP, *K, *F, *NU, *rho, *EPS, *Tmelt, *GammaT, *LoverT;
	double *Qheat,*Qimp;
	double **CP_grid, **K1_grid, **K0_grid, **NU_grid, **EPS_grid, **KAPPA_grid, **K1perp_grid, **K0perp_grid;
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
	int include_latent_heat, include_convection;
	int running_from_command_line;
	int hardwireQ;
	double Qrho;
	int use_piecewise;
	int force_precalc;
	int use_my_envelope;
	int include_sph;
	FILE *fp,*fp2;
	double Qinner;
	double energy_deposited_outer;
	double energy_deposited_inner;
	double energy_slope;
	double outburst_duration;
	int instant_heat;
	double surfF, surfy;
	double angle_mu;
	int gpe;
	int force_cooling_bc, extra_heating;
	double extra_Q,extra_y;
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
	
	long idum=-32094034;

	// initialize EOS routines
  	EOS.init(1); 
  	EOS.X[1] = 1.0;   // we only have one species at each depth;

	// ----------------------------------------------------------------------
 	//   Set parameters

	// first, some defaults	
	mass=1.62;
	G.radius=11.2;
	int ngrid=100;
	G.include_latent_heat=0;
	G.include_convection=0;
	G.nuflag = 1;
	G.energy_slope=0.0;
	G.force_precalc=0;
	G.Qinner=-1.0;
	G.running_from_command_line=1;	
	G.outburst_duration = (1.0/24.0) * 1.0/(365.0);  // rapid heating for magnetar case	(1 hour)
	EOS.accr = 0;   // set crust composition
	EOS.use_potek_eos = 0;
	EOS.use_potek_cond = 1;
	EOS.B=0; //2.2e14;   // magnetic field in the crust   (set B>0 for magnetar case)
	EOS.gap = 1;    // 0 = no gap, normal neutrons
	EOS.kncrit = 0.0;  // neutrons are normal for kn<kncrit (to use this set gap=4)
	G.mdot=0.1;
	G.energy_deposited_outer=1.0;
	G.energy_deposited_inner=-1.0;
	G.rhob=1e15; 
	G.rhot=1e3;
	G.use_piecewise=0;
	G.Qrho=1e12;
	G.instant_heat = 0;
	G.include_sph=1;
	G.angle_mu=-1.0;
	G.gpe=0;
	G.force_cooling_bc=0;
	G.extra_heating=0;
	G.use_my_envelope=0;
	G.yt=1e12;
	G.extra_Q=0.0;
	G.extra_y=0.0;
	// do we want output or not?
	G.output=1;
	
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
	int commented=0;
	while (!feof(fp)) {   // we read the file line by line
		char *e=fgets(s1,200,fp);		
		// ignoring lines that begin with \n (blank) or with # (comments)
		// or with $ (temperature profile)
		if (!strncmp(s1,"##",2)) commented = 1-commented;
		if (strncmp(s1,"#",1) && strncmp(s1,"\n",1) && strncmp(s1,">",1) && commented==0) {
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
			if (!strncmp(s,"output",6)) G.output=x;
			if (!strncmp(s,"timetorun",9)) G.time_to_run=24.0*3600.0*x;			
			if (!strncmp(s,"toutburst",9)) G.outburst_duration=x;
			if (!strncmp(s,"piecewise",9)) G.use_piecewise=(int) x;
			if (!strncmp(s,"neutrinos",9)) G.nuflag=(int) x;
			if (!strncmp(s,"accreted",8)) EOS.accr=(int) x;
			if (!strncmp(s,"angle_mu",8)) G.angle_mu=x;
			if (!strncmp(s,"latent_heat",11)) G.include_latent_heat=(int) x;
			if (!strncmp(s,"convection",10)) G.include_convection=(int) x;			
			if (!strncmp(s,"cooling_bc",10)) G.force_cooling_bc=(int) x;
			if (!strncmp(s,"extra_heating",13)) G.extra_heating=(int) x;
			if (!strncmp(s,"energy_slope",12)) G.energy_slope=x;
			if (!strncmp(s,"potek_eos",9)) EOS.use_potek_eos=(int) x;
			if (!strncmp(s,"envelope",8)) G.use_my_envelope=(int) x;
			if (!strncmp(s,"extra_Q",7)) G.extra_Q=x;
			if (!strncmp(s,"extra_y",7)) G.extra_y=x;
		}
	}

	fclose(fp);

	if (G.yt < 10.0) G.yt=pow(10.0,G.yt);

	if (G.Qinner == -1.0) G.Qinner=EOS.Q;
	if (G.energy_deposited_inner == -1.0) G.energy_deposited_inner = G.energy_deposited_outer;
	
	if (EOS.Q>=0.0) {   	// the Q values are assigned directly in 'calculate_vars'
		G.hardwireQ=1;
		printf("Using supplied Qimp values and HZ composition and heating.\n");
	} else {
		G.hardwireQ=0;
		printf("Using Qimp, composition, and heating from the crust model.\n");
	}

	// Include dipole angular dependence in B
	// The B provided is the polar magnetic field strength
	if (G.angle_mu >= 0.0) EOS.B*=sqrt(0.75*G.angle_mu*G.angle_mu+0.25);
	printf("Magnetic field set to B=%lg\n", EOS.B);
		
		{
			char fname[40]="data/";
			if (argc >1) {
				strcat(fname,argv[1]);
				read_in_data(fname);
			} else {
//				read_in_data("data/XTEJ");
//				read_in_data("data/1731");
				read_in_data("data/1659");
			}	
		}	
//		read_in_data("data/1731");  // READ IN observed lightcurve
	//	read_in_data("data/1659");  // READ IN observed lightcurve
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
  
	if (G.output) {
  		G.fp=fopen("gon_out/out","w");
  		G.fp2=fopen("gon_out/prof","w");
  		fprintf(G.fp,"%d %lg\n",G.N+1,G.g);
	}

	//G.output = 0;
 	// set up the initial temperature profile
  	if (G.use_piecewise) set_up_initial_temperature_profile_piecewise(fname);
	else set_up_initial_temperature_profile();
	//G.output = 1;

	// calculate the cooling curve
	calculate_cooling_curve();   // not outputting the heating stage,so start at t=0

	// get chisq
	double chisq = calculate_chisq();

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
  	double timesofar=0.0,last_time_output=0.0;
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
		printf("Starting output\n");
		// output 
		if (G.output) {
			start_timing(&timer);
			for (int j=1; j<=ODE.kount; j++) output_result_for_step(j,G.fp,G.fp2,timesofar,&last_time_output);
			fflush(G.fp); fflush(G.fp2);
			stop_timing(&timer,"output");
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
	printf("chisq = %lg\n", chisq);
	printf("chisq_nu = %lg/(%d-3) = %lg\n", chisq, G.obs_n, chisq/(G.obs_n-3));
	return chisq;

}



void read_in_data(const char *fname) 
{
	if (1) {   
	
	/*	
	// hardcode the data for 1731
	double t0=51930.5;
	G.obs_n=8;
	G.obs_time = vector(1,G.obs_n);
	G.obs_temp = vector(1,G.obs_n);
	G.obs_err = vector(1,G.obs_n);
	
	G.obs_time[1]=51995.1; G.obs_temp[1]=103.2; G.obs_err[1]=1.7;
	G.obs_time[2]=52165.7; G.obs_temp[2]=88.9; G.obs_err[2]=1.3;
	G.obs_time[3]=52681.6; G.obs_temp[3]=75.5; G.obs_err[3]=2.2;
	G.obs_time[4]=52859.5; G.obs_temp[4]=73.3; G.obs_err[4]=2.3;
	G.obs_time[5]=53430.5; G.obs_temp[5]=71.0; G.obs_err[5]=1.8;
	G.obs_time[6]=53500.4; G.obs_temp[6]=66.0; G.obs_err[6]=4.5;
	G.obs_time[7]=53525.4; G.obs_temp[7]=70.3; G.obs_err[7]=2.1;
	G.obs_time[8]=54969.7; G.obs_temp[8]=63.1; G.obs_err[8]=2.1;
	
	for (int i=1; i<=G.obs_n; i++) {
		G.obs_time[i]-=t0;
	}
	*/
	
	// hardcode the data for 1659
	double t0=52159.5;
	G.obs_n=7;
	G.obs_time = vector(1,G.obs_n);
	G.obs_temp = vector(1,G.obs_n);
	G.obs_err = vector(1,G.obs_n);
	
	G.obs_time[1]=52197.8; G.obs_temp[1]= 121; G.obs_err[1]= 1;
	G.obs_time[2]=52563.2; G.obs_temp[2]= 85; G.obs_err[2]= 1;
	G.obs_time[3]=52712.2; G.obs_temp[3]= 77; G.obs_err[3]= 1;
	G.obs_time[4]=52768.9; G.obs_temp[4]= 73; G.obs_err[4]= 1;
	G.obs_time[5]=53560.0; G.obs_temp[5]= 58; G.obs_err[5]= 2;
	G.obs_time[6]=53576.7; G.obs_temp[6]= 54; G.obs_err[6]= 3;
	G.obs_time[7]=54583.8; G.obs_temp[7]= 56; G.obs_err[7]= 2;
	G.obs_time[8]=56113; G.obs_temp[8]= 48.8; G.obs_err[8]= 1.6;
	
	for (int i=1; i<=G.obs_n; i++) {
		G.obs_time[i]-=t0;
	}
	
	
} else {	
	
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
}



void start_timing(clock_t *time)
{
	*time=clock();
}

void stop_timing(clock_t *time, const char*string)
{
  	printf("Time taken for %s =%lg s\n", string, (double) (clock()-*time)/((double) CLOCKS_PER_SEC)); 	
}


void output_result_for_step(int j, FILE *fp, FILE *fp2,double timesofar,double *last_time_output) 
{
	for (int i=1; i<=1; i++) calculate_vars(i,ODE.get_y(i,j),G.P[i],&G.CP[i],&G.K[i],&G.NU[i],&G.EPS[i]);
	//	for (int i=1; i<=G.N+1; i++) calculate_vars(i,ODE.get_y(i,j),G.P[i],&G.CP[i],&G.K[i],&G.NU[i],&G.EPS[i]);
	double T0;
	outer_boundary(ODE.get_y(1,j),G.K[1],G.CP[1],G.NU[1],G.EPS[1],&T0,&G.K[0],&G.CP[0],&G.NU[0],&G.EPS[0]);
	//for (int i=1; i<=G.N+1; i++) G.F[i] = calculate_heat_flux(i,ODE.get_y(i,j));
	//G.F[i]=0.5*G.g*(G.K[i]+G.K[i-1])*(ODE.get_y(i,j)-ODE.get_y(i-1,j))/G.dx;
	//G.F[1]=0.5*(G.K[0]+G.K[1])*(ODE.get_y(1,j)-T0)/G.dx;
	double dt;
	if (j==1) dt=ODE.get_x(j); else dt=ODE.get_x(j)-ODE.get_x(j-1);


	double *TT;
	TT=vector(1,G.N+1);
	for (int i=1; i<=1; i++) TT[i]=ODE.get_y(i,j);
//	for (int i=1; i<=G.N+1; i++) TT[i]=ODE.get_y(i,j);
	double FF = calculate_heat_flux(1,TT);
	for (int i=1; i<=2; i++) G.F[i] = calculate_heat_flux(i,TT);
//	for (int i=1; i<=G.N+1; i++) G.F[i] = calculate_heat_flux(i,TT);

	//FF=G.F[1];

	//G.F[2] = calculate_heat_flux(2,TT);

	free_vector(TT,1,G.N+1);

	double Lnu=0.0;
	for (int i=1; i<=G.N; i++) Lnu += G.NU[i]*G.dx*G.P[i]/G.g;
	
	if ((fabs(log10(fabs(timesofar+ODE.get_x(j))*G.ZZ)-log10(fabs(*last_time_output))) >= 0.01) || (fabs(timesofar)+ODE.get_x(j))*G.ZZ < 1e5) {
	
	// output time, fluxes and TEFF that are already redshifted
	fprintf(fp2, "%lg %lg %lg %lg %lg %lg %lg %lg %lg %lg\n", (timesofar+ODE.get_x(j))*G.ZZ, 
			pow((G.radius/11.2),2.0)*G.F[2]/(G.ZZ*G.ZZ), pow((G.radius/11.2),2.0)*FF/(G.ZZ*G.ZZ),
			ODE.get_y(G.N-5,j), pow((G.g/2.28e14)*TEFF.get(ODE.get_y(1,j))/5.67e-5,0.25)/G.ZZ, 
			ODE.get_y(1,j), pow((G.g/2.28e14)*TEFF.get(ODE.get_y(1,j))/5.67e-5,0.25),
			pow((G.radius/11.2),2.0)*G.F[G.N+1]/(G.ZZ*G.ZZ),pow((G.radius/11.2),2.0)*G.F[G.N]/(G.ZZ*G.ZZ),
			4.0*PI*pow(1e5*G.radius,2.0)*Lnu/(G.ZZ*G.ZZ));

		if ((fabs(log10(fabs(timesofar+ODE.get_x(j))*G.ZZ)-log10(fabs(*last_time_output))) >= 1000.0) || (fabs(timesofar)+ODE.get_x(j))*G.ZZ < 1e9) {
//			if (j % 1 == 0 || j==ODE.kount) {   // output every nth cycle
		// temperature profile
		fprintf(fp,"%lg\n",G.ZZ*(timesofar+ODE.get_x(j)));
		double del,gamma;
		for (int i=1; i<=G.N+1; i++) {      
			EOS.P=G.P[i]; EOS.T8=1e-8*ODE.get_y(i,j); EOS.rho=G.rho[i];//EOS.find_rho();
			if (i>1) del = (ODE.get_y(i+1,j)-ODE.get_y(i-1,j))/(2.0*G.dx*ODE.get_y(i,j));
			else del=1.0;
			gamma = G.GammaT[i]/(1e8*EOS.T8);
			EOS.set_comp();
			//if (G.ZZ*(timesofar+ODE.get_x(j)) > 0.0)
				fprintf(fp, "%lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg\n", 
				G.P[i], ODE.get_y(i,j), G.F[i], G.NU[i],
				G.g*(G.F[i+1]-G.F[i])/(G.dx*G.P[i]), EOS.rho, EOS.CP()*EOS.rho, 
				ODE.get_d(i,j),1e8*pow(EOS.P/2.521967e17,0.25), EOS.opac(), del, EOS.del_ad(),
				2.521967e-15*pow(ODE.get_y(i,j),4)/G.P[i], gamma,EOS.eps_nu());
					//G.ZZ*(timesofar+ODE.get_x(j)));	
 		}
	}
		*last_time_output=(timesofar+ODE.get_x(j))*G.ZZ;
	}
}




// --------------------------------------------------------------------------
// derivatives for the time-dependent code

void derivs(double t, double T[], double dTdt[])
// calculates the time derivatives for the whole grid
{
	// First calculate quantities at each grid point
	for (int j=1; j<=G.N; j++) calculate_vars(j,T[j], G.P[j], &G.CP[j], &G.K[j], &G.NU[j],&G.EPS[j]);
  	outer_boundary(T[1],G.K[1],G.CP[1],G.NU[1],G.EPS[1],&T[0],&G.K[0],&G.CP[0],&G.NU[0],&G.EPS[0]);
  	inner_boundary(T[G.N],G.K[G.N],G.CP[G.N],G.NU[G.N],G.EPS[G.N],
				&T[G.N+1],&G.K[G.N+1],&G.CP[G.N+1],&G.NU[G.N+1],&G.EPS[G.N+1]);

	// determine the fluxes at the half-grid points
	//  G.F[i] is the flux at i-1/2
  	for (int i=1; i<=G.N+1; i++)   G.F[i] = calculate_heat_flux(i,T);	

	// find the ocean-crust boundary: I've tried two different ways, either moving in
	// or out -- it's an issue when there is a jump in Tmelt in the ocean, e.g.
	// because of an electron capture
	int imelt;
	if (0) {
		// (i) imelt is set to be the first zone going outwards where Gamma at i-1/2 is <175
		imelt=G.N;
		while (imelt--, G.GammaT[imelt-1]/(0.5*(T[imelt]+T[imelt-1])) >= 2*175.0);
	} else {
		// (ii) imelt is set to be the first zone going inwards where Gamma at i+1/2 is >175
		imelt=1;
		while (imelt++, G.GammaT[imelt]/(0.5*(T[imelt]+T[imelt+1])) <= 2*175.0);
	}
	
	// include the latent heat as an additional term in the heat capacity at the boundary
	if (G.include_latent_heat) {
		double Gp = G.GammaT[imelt]/(0.5*(T[imelt]+T[imelt+1]));
		double Gm = G.GammaT[imelt-1]/(0.5*(T[imelt]+T[imelt-1]));
		G.CP[imelt]+= 4.0*G.LoverT[imelt] * pow(175.0/Gm,4.0)/ (pow(Gp/Gm,4.0)-1.0);
	}
	
	// include convective fluxes (only if we are cooling)
	if (G.include_convection && G.P[imelt]/G.g>1e9) {

		if (!G.accreting) {
			
			// This is the simple F=Ay prescription for the convective flux from Zach
			double AA=2.4e9;
			// first try -> double AA=2.4e11;

			AA=1e9;

			// First calculate dT/dt for the zone containing the liquid/solid boundary,
			// we need this to get d y_L/ dt which determines the convective flux
			G.CP[imelt] += AA/G.dx;
			dTdt[imelt] = G.g*(G.F[imelt+1]-G.F[imelt])/(G.dx*G.CP[imelt]*G.P[imelt]);

			// If the boundary is moving in while cooling, then there shouldn't be any convection, 
			// so turn it off. Although this doesn't seem to affect the lightcurves..
			if (dTdt[imelt] > 0.0) {
				G.CP[imelt] -= AA/G.dx;
				dTdt[imelt] = G.g*(G.F[imelt+1]-G.F[imelt])/(G.dx*G.CP[imelt]*G.P[imelt]);
				AA=0.0;
			}

			// Add in the convective flux where liquid
			for (int i=2; i<imelt; i++) {
				double P2 = exp(log(G.P[i])+0.5*G.dx);
				if (P2 > 1e25) G.F[i+1]+=AA*P2*dTdt[imelt]/G.g;	
					//else G.F[i+1]+=-G.mdot*8.8e4*9.64e17*0.02*pow(G.rho[i]*1e-9,0.333);
			}
		} else {
			
			// Add in the convective flux where liquid
			for (int i=2; i<imelt; i++) {
				double P2 = exp(log(G.P[i])+0.5*G.dx);
				if (P2 > 1e20) G.F[i+1]+=-G.mdot*8.8e4*9.64e17*0.02*pow(G.rho[i]*1e-9,0.333);
			}
		}

	}
	
	// Calculate the derivatives dT/dt
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
	if (i>1 || (G.accreting && G.outburst_duration > 1.0/365.0 && !G.force_cooling_bc))
//		if (i>1 || (G.accreting && EOS.B == 0.0))   
		// use this inside the grid, or at the surface when we are accreting (which 
		// fixes the outer temperature)
	 	flux = 0.5*(G.K[i]+G.K[i-1])*(T[i]-T[i-1])/G.dx;	
	else {
		// cooling boundary condition
		if (EOS.B == 0.0 || G.use_my_envelope) {
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
				//chi2 = sqrt(1.0 + 0.1076*B12*pow(0.03+T9,-0.559))/
				//			pow(1.0+0.819*B12/(0.03+T9),0.6463);
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
// calculates the time derivative for grid cell i 
// This is used when calculating the jacobian in tri-diagonal form
// --- i.e. as long as include_convection is not set
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

	return f;
}

void outer_boundary(double T1, double K1, double CP1, double NU1, double EPS1,
		double *T0, double *K0, double *CP0, double *NU0, double *EPS0)  
{
	if (G.accreting && G.Tt>0.0) *T0=G.Tt;   // constant temperature during accretion
	else *T0=T1*(8.0-G.dx)/(8.0+G.dx);   // assumes radiative zero solution, F\propto T^4
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
  	double e=0.01;

	if (G.include_convection) {
	// This version calculates the full set of derivatives each time
		double *f;
		f=vector(1,n);
  		for (int i=1; i<n; i++) {
    		if (i>1) {
				T[i-1]*=1.0+e; derivs(t,T,f);
      			T[i-1]/=1.0+e; dfdT[i][i-1]=(f[i]-dfdt[i])/(T[i-1]*e);
    		}
			T[i]*=1.0+e; derivs(t,T,f);
    		T[i]/=1.0+e; dfdT[i][i]=(f[i]-dfdt[i])/(T[i]*e);
    		if (i<=n) {
				T[i+1]*=1.0+e; derivs(t,T,f);
      			T[i+1]/=1.0+e; dfdT[i][i+1]=(f[i]-dfdt[i])/(T[i+1]*e);
    		}
  		}
		free_vector(f,1,n);
		
	} else {
		// takes advantage of the tri-diagonal nature to calculate as few dTdt's as needed
		double f;
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
}  


void calculate_vars(int i, double T, double P, double *CP, double *K, double *NU,double *EPS)
{
	// sometimes we get a nan value for T here from the integrator
	// In this case, set the temperature to be some value.. this seems to
	// deal with this problem ok
	if (isnan(T) || T<0.0) T=1e7;
	
	double beta=log10(T);
	
	if (beta > G.betamax) beta = G.betamax;
	if (beta < G.betamin) beta = G.betamin;
		
	if (beta <= G.betamax && beta >= G.betamin) {
		// lookup values in the precalculated table
		int j = 1 + (int) ((beta-G.betamin)/G.deltabeta);
		double interpfac=(beta-(G.betamin + (j-1)*G.deltabeta))/G.deltabeta;
		// interpolate the thermal conductivity to the current
		// value of impurity parameter Q
		double K0,K1,K0perp,K1perp;
		K0=G.K0_grid[i][j] + (G.K0_grid[i][j+1]-G.K0_grid[i][j])*interpfac;
		K1=G.K1_grid[i][j] + (G.K1_grid[i][j+1]-G.K1_grid[i][j])*interpfac;
		K0perp=G.K0perp_grid[i][j] + (G.K0perp_grid[i][j+1]-G.K0perp_grid[i][j])*interpfac;
		K1perp=G.K1perp_grid[i][j] + (G.K1perp_grid[i][j+1]-G.K1perp_grid[i][j])*interpfac;
		//K0perp=0.0; K1perp=0.0;

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

		double kappa;
		kappa=G.KAPPA_grid[i][j] + (G.KAPPA_grid[i][j+1]-G.KAPPA_grid[i][j])*interpfac;
		kappa*=G.g;
		KK += kappa;
		
		if (EOS.B > 0) {
			KKperp=0.0; //G.g*K0perp*K1perp/(K0perp*Qval+(1.0-Qval)*K1perp);
			if (G.angle_mu >= 0.0) {
				KK *= 4.0*G.angle_mu*G.angle_mu/(1.0+3.0*G.angle_mu*G.angle_mu);
			} else {
				KK = 0.5*(1.0544*KK+0.9456*KKperp);  // average over dipole geometry	
			}
		}
//		if (EOS.B > 0) {
//		//KKperp = G.g*K0perp*K1perp/(K0perp*Qval+(1.0-Qval)*K1perp);		
//		double fcond = 4.0*G.angle_mu*G.angle_mu/(1.0+3.0*G.angle_mu*G.angle_mu);		
//		*K=fcond*KK;//+(1.0-fcond)*KKperp;	
//	} else {
		*K=KK;
//	}
		
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
		printf("Note: log10T outside range (T=%lg beta=%lg)\n", T,beta);
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
			//Kcondperp=0.0;
			if (EOS.B > 0.0) {
			if (G.angle_mu >= 0.0) {
				Kcond *= 4.0*G.angle_mu*G.angle_mu/(1.0+3.0*G.angle_mu*G.angle_mu);
			} else {
				Kcond = 0.5*(1.0544*Kcond+0.9456*Kcondperp);  // average over dipole geometry
			}
			}
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
	int commented=0;
	while (!feof(fp)) {		
		double rho, T;
		// new lines: read from lines marked ">" in init.dat
		char *e=fgets(s1,200,fp);		
		if (!strncmp(s1,"##",2)) commented = 1-commented;
		if (!strncmp(s1,">",1) && commented==0) {
			sscanf(s1,">%lg\t%lg\n",&rho,&T);
			// old: direct read from Tinit.dat
			//		fscanf(fp,"%lg %lg\n",&rho,&T);
			if (rho>=-1.0) {
				if (rho == 0.0) rho=G.rho[1];
				if (rho == -1.0) rho=G.rho[G.N];
				rhovec[i] = rho;
				if (T < 0.0) T=G.Tc;
				Tvec[i] = T;
				//printf("%lg %lg\n",rhovec[i],Tvec[i]);
				i++;
			}
		}
	}
	fclose(fp);
	int nvec=i-1;

	if (nvec == 0) {
		printf("ERROR:The piecewise flag is set but the temperature profile is not specified in init.dat!\n");
		exit(0);
	}

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
				while (rhovec[j] < G.rho[i] && j<nvec) j++;
				Ti = pow(10.0,log10(Tvec[j-1]) + log10(Tvec[j]/Tvec[j-1])*log10(G.rho[i]/rhovec[j-1])/log10(rhovec[j]/rhovec[j-1]));
				//printf("%d %lg %lg %lg %lg %lg\n",j,Tvec[j-1],Tvec[j],rhovec[j-1],rhovec[j],G.rho[i]);
			}
		}	
		
//		printf("%lg %lg\n",G.rho[i],Ti);
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
		//if (Ed > 0.0) printf("heating cell %d: Ti=%lg Tf=%lg E25=%lg rho=%lg\n",
		//			i, G.Tc, Ti, Ed*G.rho[i]*1e-25, G.rho[i]);
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
			
			fprintf(fp, "%d %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg\n", i, G.P[i], Ti, EOS.rho,EOS.CV(), 
					Kcond, EOS.Yn, 1e-39*EOS.Yn * EOS.rho/1.67e-24, tt,E,EOS.A[1],EOS.Z[1],TTC, EOS.econd(), EOS.Ye()*EOS.rho/1.67e-24);
		}
	}	
	
	printf("Total energy input to get this initial T profile = %lg (redshifted=%lg)\n",totalEd,totalEd/G.ZZ);
	
	if (G.output) fclose(fp);
	
	free_vector(rhovec,1,100);
	free_vector(Tvec,1,100);
	
}





void set_up_initial_temperature_profile(void)
// initialize the temperature profile on the grid
{
	
	// first run some accretion with heating, long enough for the crust
	// to get into a thermal steady-state
	for (int i=G.N+1; i>=1; i--) {
		// a linear profile between top and bottom
		//double Ti = pow(10.0,log10(G.Tc) + log10(0.3*G.Tt/G.Tc)*log10(G.P[i]/G.Pb)/log10(G.Pt/G.Pb));
		// or constant profile
		double Ti = G.Tc;
		// a linear profile adjusts to steady state *much* more quickly,
		// but for XTEJ for example I want to heat up from isothermal and the crust 
		// does not get to steady state
		ODE.set_bc(i,Ti);
	}
	clock_t timer;
	start_timing(&timer);

	if (1) {
	// First, let the crust cool for 30 years to get into eqm with the core
	G.accreting = 0;  // switch off heating for this one
	ODE.go(0.0, 30.0*3.15e7, 1e6, 1e-6, derivs);
	for (int i=1; i<=G.N+1; i++)
		ODE.set_bc(i,ODE.get_y(i,ODE.kount));
	}

	G.accreting = !G.instant_heat;	
	double dt=G.outburst_duration*3.15e7*0.0001;
	if (dt > 1e6) dt=1e6;
	//	ODE.verbose  = 1;
//	ODE.dxsav = 1e8;
	ODE.go(0.0, G.outburst_duration*3.15e7,dt, 1e-6, derivs);
	stop_timing(&timer,"ODE.go (initial heating)");
	printf("number of steps = %d  (time=%lg)\n", ODE.kount, ODE.get_x(ODE.kount));

	double timesofar = -G.outburst_duration*3.15e7;
	double last_time_output = timesofar;
	if (G.output) {
		for (int j=1; j<=ODE.kount; j++) output_result_for_step(j,G.fp,G.fp2,timesofar,&last_time_output);
		fflush(G.fp); fflush(G.fp2);
	}

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
			printf(" Tf=%lg  rho=%lg\n", Ti, G.rho[i]);
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
			//double Kcond = potek_cond(&Kcondperp);
			double Kcond = EOS.K_cond(EOS.Chabrier_EF());
								
			I+=sqrt(EOS.CV()/(Kcond*EOS.rho))*G.P[i]*G.dx/G.g;
			double tt = I*I*0.25/(24.0*3600.0);
	//		tt = 0.25 * EOS.rho * EOS.CP() * pow(G.P[i]/(G.g*EOS.rho),2.0)/(Kcond*24.0*3600.0);

			E+=energy_deposited(i)*crust_heating(i)*G.mdot*G.g*G.outburst_duration*3.15e7*
				4.0*PI*G.radius*G.radius*1e10*G.dx*G.P[i]/G.g;
											
			fprintf(fp, "%d %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg\n", i, G.P[i], Ti, EOS.rho,EOS.CV(), 
					Kcond, EOS.Yn, 1e-39*EOS.Yn * EOS.rho/1.67e-24, tt,E,EOS.A[1],EOS.Z[1],EOS.TC(),
					EOS.Chabrier_EF(),Kcondperp,0.4*1e-9*EOS.rho*pow(EOS.Ye()/0.4,3.0)*EOS.Z[1]/34.0, EOS.cve, EOS.cvion,
					EOS.chi(&EOS.rho), EOS.chi(&EOS.T8), G.P[i]/(G.g*EOS.rho), EOS.econd());
				EOS.Q=Q_store;	
		}
	}	
	if (G.output) fclose(fp);

	printf("Total energy deposited=%lg (redshifted=%lg)\n",E,E/G.ZZ);

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
	G.KAPPA_grid = matrix(1,G.N+2,1,G.nbeta);	
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
				//Kcond = 3.03e20*pow(EOS.T8,3)/(EOS.opac()*G.y[i]);
				//Kcondperp=Kcond;
				//if (Kcondperp < 1e-2*Kcond) Kcondperp=1e-2*Kcond;
				G.K0_grid[i][j]=EOS.rho*Kcond/G.P[i];
				G.K0perp_grid[i][j]=EOS.rho*Kcondperp/G.P[i];
				EOS.Q=1.0;
				//Kcond = EOS.K_cond(EOS.Chabrier_EF());
				//Kcond = 3.03e20*pow(EOS.T8,3)/(EOS.opac()*G.y[i]);
				//Kcondperp=Kcond;
				Kcond = potek_cond(&Kcondperp);
				//if (Kcondperp < 1e-2*Kcond) Kcondperp=1e-2*Kcond;
				G.K1_grid[i][j]=EOS.rho*Kcond/G.P[i];
				G.K1perp_grid[i][j]=EOS.rho*Kcondperp/G.P[i];
				EOS.Q=Q_store;  // restore to previous value

				double kappa=EOS.opac();
				G.KAPPA_grid[i][j] = 3.03e20*pow(EOS.T8,3)/(EOS.kappa_rad*G.P[i]);
				//G.KAPPA_grid[i][j] = 0.0;

				// We include thermal conductivity only. The equivalent expression but putting in all sources
				// of opacity would be
				//3.03e20*pow(EOS.T8,3)/(EOS.opac()*G.y[i]);

				fprintf(fp, "%lg %lg %lg %lg %lg %lg %lg %lg %lg\n", EOS.T8, G.CP_grid[i][j], 
					G.K0_grid[i][j],G.K1_grid[i][j], G.K0perp_grid[i][j],G.K1perp_grid[i][j],
					G.NU_grid[i][j], G.EPS_grid[i][j], G.KAPPA_grid[i][j] );

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
				fscanf(fp, "%lg %lg %lg %lg %lg %lg %lg %lg %lg\n", &EOS.T8, &G.CP_grid[i][j], 
					&G.K0_grid[i][j],&G.K1_grid[i][j], &G.K0perp_grid[i][j],&G.K1perp_grid[i][j],
					&G.NU_grid[i][j], &G.EPS_grid[i][j],&G.KAPPA_grid[i][j]);
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
	double ener;
	if (G.rho[i]>4e11) ener = G.energy_deposited_inner;
	else ener = G.energy_deposited_outer;
	ener *= pow(G.rho[i]/1e10,G.energy_slope);
	return ener;
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

		// Put in a power-law heating:
//		if (G.rho[i]>1e9) 
//		eps *= pow(G.rho[i]/1e9,1.0);
//		eps *= G.rho[i]/1e11;
		
		// the next line limits the heating to a limited region of the crust
		if (G.rho[i] < G.rhot) eps=0.0;
		if (G.rho[i] > G.rhob) eps=0.0;
		
	} else {  // otherwise we are doing an accreting neutron star

		if (!G.hardwireQ) 
			eps = G.Qheat[i]*8.8e4*9.64e17/(G.P[i]*G.dx);
		else
			eps = crust_heating_rate(i);
	}

	return eps;	
}


double crust_heating_rate(int i) 
// calculates the crust heating in erg/g/s
 // except for a factor of gravity --- multiply by gravity to get these units
{
	double eps=0.0, P = G.P[i];
	// simple "smeared out" heating function, 1.2MeV in inner crust, 0.2MeV in outer crust
	if (P >= 1e16*2.28e14 && P <= 1e17*2.28e14) eps=8.8e4*1.7*9.64e17/(P*log(1e17/1e16));
// 	if (y >= 1e12 && y < 1e15) eps=G.mdot*8.8e4*0.15*9.64e17/(y*log(1e15/1e12));
 	if (P >= 3e12*2.28e14 && P < 3e15*2.28e14) eps=8.8e4*0.2*9.64e17/(P*log(3e15/3e12));
//if (y >= 6e15 && y <= 3e18) eps=G.mdot*8.8e4*1.2*9.64e17/(y*log(3e18/6e15));

	// Extra heat source in the ocean
	if (G.extra_heating) {	
	
		if (G.P[i]*exp(-0.5*G.dx) <G.extra_y*2.28e14 && G.P[i]*exp(0.5*G.dx)>G.extra_y*2.28e14)
				eps+=8.8e4*G.extra_Q*9.64e17/(P*G.dx);
	
//	 && P >=0.9*G.extra_y*2.28e14 && P<=2.0*G.extra_y*2.28e14) eps+=8.8e4*G.extra_Q*9.64e17/(P*G.dx);
//	if (G.extra_heating && P >=0.5*G.extra_y*2.28e14 && P<=2.0*G.extra_y*2.28e14) eps+=8.8e4*G.extra_Q*9.64e17/(P*log(3.0/0.333));
	}
	
	
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
//	if (G.use_my_envelope) fp = fopen("out/grid","r");
	if (G.use_my_envelope) {
		if (EOS.B == 1e15) fp = fopen("out/grid_1e15_nopotek","r");
		else if (EOS.B == 1e14) fp = fopen("out/grid_1e14_potek","r");
		else if (EOS.B == 3e14) fp = fopen("out/grid_3e14_potek","r");
		else if (EOS.B == 3e15) fp = fopen("out/grid_3e15_potek","r");
		else {
			printf("Don't know which envelope model to use for this B!\n");
			exit(1);
		}
	}
	else {
		//fp = fopen("out/grid","r");
		if (G.gpe) fp = fopen("out/grid_He4","r");
		else fp = fopen("out/grid_He9","r");
	}
	//FILE *fp = fopen("grid_sorty","r");
	FILE *fp2 = fopen("gon_out/TbTeff", "w");
	
	double y,T,F,rho,dummy;
	int count = 0;
	while (!feof(fp)) {
		
		fscanf(fp, "%lg %lg %lg %lg %lg %lg\n", &y, &T, &F,&rho,&dummy,&dummy);
		//printf("%lg %lg %lg\n", y, T, F);
		if (fabs(y-log10(G.yt))<1e-3) {  // select out the points which correspond to the top column
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
//  	G.Pt=G.yt*G.g;   // pressure at the top
	G.Pt=G.yt*2.28e14;   // pressure at the top   // note I need to use 2.28 here to get the correct match to the envelope

	Spline QiSpline;
	Spline QhSpline;
	if (!G.hardwireQ) {   // only need this if we read in the crust model
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
		dens[i] = 5.0 + (i-1)*(15.0-5.0)*0.001;
		EOS.rho = pow(10.0,dens[i]);
		set_composition();
		pres[i] = log10(EOS.ptot());
		//printf("%d %lg %lg %lg %lg %lg\n",i,EOS.rho,dens[i],pres[i],EOS.A[1],EOS.Z[1],EOS.Yn);
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

	double Qtot=0.0;
  	for (int i=0; i<=G.N+2; i++) {
    	double x=log(G.Pt)+G.dx*(i-1);
    	G.P[i]=exp(x);

		double P1 = exp(x-0.5*G.dx);
		double P2 = exp(x+0.5*G.dx);
		double rho1 = pow(10.0,RHO.get(log10(P1)));
		double rho2 = pow(10.0,RHO.get(log10(P2)));
		G.Qheat[i]=0.0;
		if (!G.hardwireQ) {
			G.Qheat[i]=QhSpline.get(log10(rho2))-QhSpline.get(log10(rho1));
			if (G.Qheat[i]<0.0) G.Qheat[i]=0.0;
		}
		Qtot+=G.Qheat[i];

		if (!G.hardwireQ) {
			G.Qimp[i]=QiSpline.get(log10(G.rho[i]));
			if (G.Qimp[i] < 0.0) G.Qimp[i]=0.0;
		}
	
		EOS.rho = rho2;
		set_composition();
		// GammaT[i] refers to i+1/2
		// The following line uses a composition of 56Fe to calculate gamma,
		// it avoids jumps in the melting point associated with e-capture boundaries
		if (0) {
			double Z=26.0,A=56.0;   // 28Si
			G.GammaT[i] = pow(Z*4.8023e-10,2.0)*pow(4.0*PI*EOS.rho/(3.0*A*1.67e-24),1.0/3.0)/1.38e-16;
		} else {
		G.GammaT[i] = pow(EOS.Z[1]*4.8023e-10,2.0)*pow(4.0*PI*EOS.rho/(3.0*EOS.A[1]*1.67e-24),1.0/3.0)/1.38e-16;
		}
		
    	G.rho[i] = pow(10.0,RHO.get(log10(G.P[i])));
		EOS.rho = G.rho[i];
		set_composition();
		
		G.Tmelt[i] = 5e8*pow(G.P[i]/(2.28e14*1.9e13),0.25)*pow(EOS.Z[1]/30.0,5.0/3.0);
		G.LoverT[i] = 0.8 * 1.38e-16 /(EOS.A[1]*1.67e-24);
		fprintf(fp, "%d %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg\n", i, G.P[i], G.rho[i], EOS.A[1]*(1.0-EOS.Yn), 
			EOS.Z[1], EOS.Yn,EOS.A[1],EOS.ptot(), G.Tmelt[i], G.GammaT[i]/1e8, G.LoverT[i]*1e8);
			//,AASpline.get(log10(G.rho[i])),  ZZSpline.get(log10(G.rho[i])), G.Qimp[i], G.Qheat[i]);
  	}

	fclose(fp);

	if (!G.hardwireQ) QiSpline.tidy();

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



void heatderivs2(double T, double E[], double dEdT[])
// calculates the time derivatives for the whole grid
{
	EOS.T8 = T/1e8;
	EOS.rho = EOS.find_rho();
	dEdT[1] = EOS.CP();
}

void heatderivs(double E, double T[], double dTdE[])
// calculates the time derivatives for the whole grid
{
	EOS.T8 = T[1]/1e8;
	EOS.rho = EOS.find_rho();
	dTdE[1] = 1e25/(EOS.rho*EOS.CP());
}

