// crustcool.cc
//

#include <stdio.h>
#include <string.h>
#include <math.h>
#include "../h/spline.h"
#include "../h/nrutil.h"
#include "../h/odeint.h"
#include "../h/crust.h"
#include "../h/data.h"

void parse_parameters(char *fname,char *sourcename,Crust *crust,double *,int*,int*);
void set_up_initial_temperature_profile_piecewise(char *fname, Crust *crust);


int main(int argc, char *argv[])
{
	// Initialize the crust
	Crust crust;
	
	// Get input parameters
	// determine the filename for the 'init.dat' parameter file
	char fname[200];
	char fnamedefault[10]="init.dat";
	switch(argc) {
		case 3:
			strcat(fname,"/tmp/init.dat.");
			strcat(fname,argv[1]);
			break;
		case 2:
			strcat(fname,"init/init.dat.");
			strcat(fname,argv[1]);
			break;
		default:
			strcat(fname,fnamedefault);
	}
	char sourcename[200]="1659";
	double time_to_run=1e4;
	int instant_heat=0, use_piecewise=0;
	parse_parameters(fname,sourcename,&crust,&time_to_run,&instant_heat,&use_piecewise);
			
	// Setup
	crust.setup();
	
	// Evolve the crust in time
	
	// Heating phase 
	if (use_piecewise) {
		set_up_initial_temperature_profile_piecewise(fname,&crust);
	} else {
		crust.output=0;
		crust.evolve(crust.outburst_duration*365.0,crust.mdot);
	}

	// Cooling phase
	crust.output=1;
	crust.force_precalc=0;
	crust.evolve(time_to_run,0.0);
	
	// Calculate the chi-sq
	read_in_data(sourcename);
	calculate_chisq(&crust.ODE,&crust.TEFF,crust.g,crust.ZZ,crust.radius,crust.Lscale,crust.Lmin);
}





void parse_parameters(char *fname,char *sourcename,Crust *crust,double *time_to_run,int *instant_heat,int*use_piecewise) {
 	// Set parameters
	printf("============================================\n");
	printf("Reading input data from %s\n",fname);
	FILE *fp = fopen(fname,"r");
	char s1[100];
	char s[100];
	double x;				
	int commented=0;
	while (!feof(fp)) {   // we read the file line by line
		(void) fgets(s1,200,fp);		
		// ignoring lines that begin with \n (blank) or with # (comments)
		// or with $ (temperature profile)
		if (!strncmp(s1,"##",2)) commented = 1-commented;
		if (strncmp(s1,"#",1) && strncmp(s1,"\n",1) && strncmp(s1,">",1) && commented==0) {
			sscanf(s1,"%s\t%lg\n",s,&x);
			if (!strncmp(s,"Bfield",6)) crust->B=x;
			if (!strncmp(s,"Tc",2)) crust->Tc=x;
			if (!strncmp(s,"Tt",2)) crust->Tt=x;
			if (!strncmp(s,"SFgap",5)) crust->gap=(int) x;
			if (!strncmp(s,"ngrid",5)) crust->N=(int) x;
			if (!strncmp(s,"kncrit",6)) crust->kncrit=x;
			if (!strncmp(s,"mdot",4)) crust->mdot=x;
			if (!strncmp(s,"mass",4)) crust->mass=x;
			if (!strncmp(s,"gpe",3)) crust->gpe=(int) x;
			if (!strncmp(s,"radius",6)) crust->radius=x;
			if (!strncmp(s,"Edep",4)) crust->energy_deposited_outer=x;
			if (!strncmp(s,"ytop",4)) crust->yt=x;
			if (!strncmp(s,"Einner",6)) crust->energy_deposited_inner=x;
			if (!strncmp(s,"Qimp",4)) crust->Qimp=x;
			if (!strncmp(s,"Qrho",4)) crust->Qrho=x;
			if (!strncmp(s,"rhob",4)) crust->rhob=x;
			if (!strncmp(s,"rhot",4)) crust->rhot=x;
			if (!strncmp(s,"precalc",7)) crust->force_precalc=(int) x;
			if (!strncmp(s,"instant",7)) *instant_heat=(int) x;
			if (!strncmp(s,"Qinner",6)) crust->Qinner=x;
			if (!strncmp(s,"output",6)) crust->output=x;
			if (!strncmp(s,"timetorun",9)) *time_to_run=x;			
			if (!strncmp(s,"toutburst",9)) crust->outburst_duration=x;
			if (!strncmp(s,"piecewise",9)) *use_piecewise=(int) x;
			if (!strncmp(s,"neutrinos",9)) crust->nuflag=(int) x;
			if (!strncmp(s,"accreted",8)) crust->accr=(int) x;
			if (!strncmp(s,"angle_mu",8)) crust->angle_mu=x;
			if (!strncmp(s,"cooling_bc",10)) crust->force_cooling_bc=(int) x;
			if (!strncmp(s,"extra_heating",13)) crust->extra_heating=(int) x;
			if (!strncmp(s,"deep_heating_factor",19)) crust->deep_heating_factor=x;
			if (!strncmp(s,"energy_slope",12)) crust->energy_slope=x;
			if (!strncmp(s,"potek_eos",9)) crust->use_potek_eos=(int) x;
			if (!strncmp(s,"envelope",8)) crust->use_my_envelope=(int) x;
			if (!strncmp(s,"extra_Q",7)) crust->extra_Q=x;
			if (!strncmp(s,"extra_y",7)) crust->extra_y=x;
			if (!strncmp(s,"Lscale",6)) crust->Lscale=x;
			if (!strncmp(s,"Lmin",4)) crust->Lmin=x;
			if (!strncmp(s,"source",6)) {
				sscanf(s1,"%s\t%s\n",s,sourcename);
			}
		}
	}

	fclose(fp);	
}




void set_up_initial_temperature_profile_piecewise(char *fname,Crust *crust)
{
	// first read in the specified temperature-density relation from the file
	FILE *fp;
	char s1[100];
	printf("Reading initial temperature profile from %s\n",fname);
	fp = fopen(fname,"r");
	double *rhovec, *Tvec;
	rhovec=vector(1,100);
	Tvec=vector(1,100);
	int commented=0;
	rhovec[1] = crust->rho[1];
	Tvec[1] = crust->Tc;
	int i=2;
	while (!feof(fp)) {		
		double rho, T,T2=0.0;
		// new lines: read from lines marked ">" in init.dat
		(void) fgets(s1,200,fp);		
		if (!strncmp(s1,"##",2)) commented = 1-commented;
		if (!strncmp(s1,">",1) && commented==0) {
			int nvar = sscanf(s1,">%lg\t%lg\t%lg\n",&rho,&T,&T2);
			// old: direct read from Tinit.dat
			//		fscanf(fp,"%lg %lg\n",&rho,&T);
			if (T < 0.0) T=crust->Tc;
			if (T2 < 0.0) T2=crust->Tc;
			if (rho < 0.0) rho = crust->rho[crust->N];
			if (rho == 0.0) {
				Tvec[1] = T;
			} else {
				rhovec[i] = rho;
				Tvec[i] = T;
				i++;
				if (nvar == 3) {
					rhovec[i] = rho*1.01;
					Tvec[i] = T2;
					i++;
				}
			}			
		}
	}
	fclose(fp);
	if (rhovec[i-1] != crust->rho[crust->N]) {  // if we didn't specify it in the file,
									// set the temperature of the base to the core temperature
		rhovec[i] = crust->rho[crust->N];
		Tvec[i] = crust->Tc;
		i++;
	}
	int nvec=i-1;

	if (nvec == 0) {
		printf("ERROR:The piecewise flag is set but the temperature profile is not specified in init.dat!\n");
		exit(0);
	}

	// now assign initial temperatures
	for (int i=1; i<=crust->N+1; i++) {
		double Ti;
		if (i==1) {
			Ti=Tvec[1]; 
			crust->Tt=Ti;
		} else {
			if (i==crust->N+1) {
				Ti=Tvec[nvec];		
			} else {
				int	j=1; 
				while (rhovec[j] < crust->rho[i] && j<nvec) j++;
				Ti = pow(10.0,log10(Tvec[j-1]) + log10(Tvec[j]/Tvec[j-1])*log10(crust->rho[i]/rhovec[j-1])/log10(rhovec[j]/rhovec[j-1]));
			}
		}	
		
		crust->T[i]=Ti;
	}
	
	free_vector(rhovec,1,100);
	free_vector(Tvec,1,100);
}
