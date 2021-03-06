// crustcool.cc
//

#include <stdio.h>
#include <string.h>
#include <math.h>
#include "../h/crust.h"
#include "../h/data.h"

void parse_parameters(char *fname,char *sourcename,Crust &crust,double &,int &, int&, int &);
void set_up_initial_temperature_profile_piecewise(char *fname, Crust &crust);


int main(int argc, char *argv[])
{
	// Initialize the crust
	Crust crust;
	
	// Get input parameters
	// determine the filename for the 'init.dat' parameter file
	char fname[200]="";
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

	// parse the input file
	char sourcename[200]="1659";
	double time_to_run=1e4;
	int use_piecewise=0, output_heating=0, output_cooling=1;
	parse_parameters(fname,sourcename,crust,time_to_run,use_piecewise,output_heating,output_cooling);

	printf("============================================\n");
			
	// Setup the crust
	crust.setup();
	
	// Evolve the crust in time
	
	// Heating phase 
	if (use_piecewise) {
		// initial temperature profile was specified in the init.dat file
		set_up_initial_temperature_profile_piecewise(fname,crust);
	} else {
		// evolve the crust while heating
		crust.output=output_heating;   // default is to not output the lightcurve while heating
		// now evolve for the set outburst_duration and mdot
		// if these were not specified in the init.dat file,
		// then they have default values of 1 hour and mdot=1.0 which is used for the 
		// magnetar case (rapid heating)
		crust.evolve(crust.outburst_duration*365.0,crust.mdot);
	}

	// Cooling phase
	if (time_to_run > 0.0) {
		crust.output=output_cooling;
		crust.evolve(time_to_run,0.0);
	}
	
	// Calculate the chi-sq
	Data data;
	data.read_in_data(sourcename);
	data.calculate_chisq(crust);
}





void parse_parameters(char *fname,char *sourcename,Crust &crust,double &time_to_run,int &use_piecewise,
			int &output_heating, int &output_cooling) {
 	// Set parameters
	printf("Reading input data from %s\n",fname);
	FILE *fp = fopen(fname,"r");
	char s[200];
	char s1[200];
	char s2[200];
	char includename[200]="";
	double x;				
	int commented=0;
	while (!feof(fp)) {   // we read the file line by line
		(void) fgets(s1,200,fp);		
		// ignoring lines that begin with \n (blank) or with # (comments)
		// or with $ (temperature profile)
		if (!strncmp(s1,"##",2)) commented = 1-commented;
		if (!strncmp(s1,"<",1)) {   // include another init file
			sscanf(s1,"%s\t%s\n",s,s2);
			strcat(includename,"init/init.dat.");
			strcat(includename,s2);			
			parse_parameters(includename,sourcename,crust,time_to_run,use_piecewise,output_heating,output_cooling);
		}
		if (strncmp(s1,"#",1) && strncmp(s1,"\n",1) && strncmp(s1,">",1) && commented==0) {
			sscanf(s1,"%s\t%lg\n",s,&x);
			if (!strncmp(s,"Bfield",6)) crust.B=x;
			if (!strncmp(s,"Tc",2)) crust.Tc=x;
			if (!strncmp(s,"Tt",2)) crust.Tt=x;
			if (!strncmp(s,"SFgap",5)) crust.gap=(int) x;
			if (!strncmp(s,"ngrid",5)) crust.N=(int) x;
			if (!strncmp(s,"kncrit",6)) crust.kncrit=x;
			if (!strncmp(s,"mdot",4)) crust.mdot=x;
			if (!strncmp(s,"mass",4)) crust.mass=x;
			if (!strncmp(s,"C_core",6)) crust.C_core=x;
			if (!strncmp(s,"Lnu_core_alpha",14)) crust.Lnu_core_alpha=x;
			if (!strncmp(s,"Lnu_core_norm",13)) crust.Lnu_core_norm=x;
			if (!strncmp(s,"gpe",3)) crust.gpe=(int) x;
			if (!strncmp(s,"radius",6)) crust.radius=x;
			if (!strncmp(s,"Edep",4)) crust.energy_deposited_outer=x;
			if (!strncmp(s,"ytop",4)) crust.yt=x;
			if (!strncmp(s,"Einner",6)) crust.energy_deposited_inner=x;
			if (!strncmp(s,"resume",6)) crust.resume=x;
			if (!strncmp(s,"Qimp",4)) crust.Qimp=x;
			if (!strncmp(s,"Qrho",4)) crust.Qrho=x;
			if (!strncmp(s,"rhob",4)) crust.rhob=x;
			if (!strncmp(s,"rhot",4)) crust.rhot=x;
			if (!strncmp(s,"precalc",7)) crust.force_precalc=(int) x;
			if (!strncmp(s,"Qinner",6)) crust.Qinner=x;
			if (!strncmp(s,"output_cooling",14)) output_cooling=x;
			if (!strncmp(s,"output_heating",14)) output_heating=x;
			if (!strncmp(s,"timetorun",9)) time_to_run=x;			
			if (!strncmp(s,"toutburst",9)) crust.outburst_duration=x;
			if (!strncmp(s,"piecewise",9)) use_piecewise=(int) x;
			if (!strncmp(s,"neutrinos",9)) crust.nuflag=(int) x;
			if (!strncmp(s,"accreted",8)) crust.accr=(int) x;
			if (!strncmp(s,"angle_mu",8)) crust.angle_mu=x;
			if (!strncmp(s,"cooling_bc",10)) crust.force_cooling_bc=(int) x;
			if (!strncmp(s,"extra_heating",13)) crust.extra_heating=(int) x;
			if (!strncmp(s,"deep_heating_factor",19)) crust.deep_heating_factor=x;
			if (!strncmp(s,"energy_slope",12)) crust.energy_slope=x;
			if (!strncmp(s,"potek_eos",9)) crust.use_potek_eos=(int) x;
			if (!strncmp(s,"envelope",8)) crust.use_my_envelope=(int) x;
			if (!strncmp(s,"extra_Q",7)) crust.extra_Q=x;
			if (!strncmp(s,"extra_y",7)) crust.extra_y=x;
			if (!strncmp(s,"Lscale",6)) crust.Lscale=x;
			if (!strncmp(s,"Lmin",4)) crust.Lmin=x;
			if (!strncmp(s,"source",6)) {
				sscanf(s1,"%s\t%s\n",s,sourcename);
			}
		}
	}

	fclose(fp);	
}




void set_up_initial_temperature_profile_piecewise(char *fname,Crust &crust)
{
	// first read in the specified temperature-density relation from the file
	FILE *fp;
	char s1[100];  // string to hold each line of the file
	printf("Reading initial temperature profile from %s\n",fname);
	fp = fopen(fname,"r");
	double rhovec[101],Tvec[101];  // temperature and density points 
	rhovec[0]=0.0; Tvec[0]=0.0;   // send "0.0" for the top density and temperature
								// Crust converts this to the right values
	int commented=0;
	int i=1;
	while (!feof(fp)) {	
		double rho, T,T2=0.0;
		// new lines: read from lines marked ">" in init.dat
		(void) fgets(s1,200,fp);
		if (!strncmp(s1,"##",2)) commented = 1-commented;
		if (!strncmp(s1,">",1) && commented==0) {
			int nvar = sscanf(s1,">%lg\t%lg\t%lg\n",&rho,&T,&T2);
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
	fclose(fp);
	int nvec=i;
	if (nvec == 1) {
		printf("ERROR:The piecewise flag is set but the temperature profile is not specified in init.dat!\n");
		exit(0);
	}
	
	crust.set_temperature_profile(rhovec,Tvec,nvec);

}
