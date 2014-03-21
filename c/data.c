#include <stdio.h>
#include <math.h>
#include "../h/spline.h"
#include "../h/nr.h"
#include "../h/nrutil.h"
#include "../h/odeint.h"
#include <string.h>

struct Data {
	double *t, *TT, *Te;
	int n;
	int luminosity;
} data;


void read_in_data(const char *sourcename) 
{
	printf("Source = %s\n",sourcename);
	data.luminosity = 0;
	
	if (1) {    // hardcoded data    

		double t0=0.0;

		if (!strncmp(sourcename,"1659",4)) {  // MXB 1659	
			t0=52159.5;
			data.n=8;    // 7 data points in BC09
			data.t = vector(1,data.n);
			data.TT = vector(1,data.n);
			data.Te = vector(1,data.n);

			data.t[1]=52197.8; data.TT[1]= 121; data.Te[1]= 1;
			data.t[2]=52563.2; data.TT[2]= 85; data.Te[2]= 1;
			data.t[3]=52712.2; data.TT[3]= 77; data.Te[3]= 1;
			data.t[4]=52768.9; data.TT[4]= 73; data.Te[4]= 1;
			data.t[5]=53560.0; data.TT[5]= 58; data.Te[5]= 2;
			data.t[6]=53576.7; data.TT[6]= 54; data.Te[6]= 3;
			data.t[7]=54583.8; data.TT[7]= 56; data.Te[7]= 2;
			data.t[8]=56113; data.TT[8]= 48.8; data.Te[8]= 1.6;
		}
		
		for (int i=1; i<=data.n; i++) {
			data.t[i]-=t0;
		}
	
	} else {	
	
		char fname[100]="data/";
		strcat(fname,sourcename);
	
		printf("Reading data from %s\n", fname);
	
		FILE *fp = fopen(fname,"r");
	
		double t0;
		fscanf(fp, "%lg %d\n", &t0, &data.n);
	
		data.t = vector(1,data.n);
		data.TT = vector(1,data.n);
		data.Te = vector(1,data.n);
		
		for (int i=1; i<=data.n; i++) {
			double dummy;
			fscanf(fp,"%lg %lg %lg %lg %lg\n",  &data.t[i], &dummy, &dummy, &data.TT[i],&data.Te[i]);
			data.t[i]-=t0;
		}	
		fclose(fp);	
	}
}



void calculate_chisq(Ode_Int *ODE, Spline *TEFF, double g, double ZZ, double R,double Lscale,double Lmin)	
// uses the result of the cooling to calculate chi-squared
{
	int nmodel = ODE->kount;
	
	// set up a spline which has the prediction for observed Teff vs time 
	Spline TE;
	double *yy = vector(1,nmodel);
	double *xx = vector(1,nmodel);
	for (int k=1; k<=nmodel; k++) { 
		xx[k]=ODE->get_x(k)*ZZ/(3600.0*24.0);
		if (data.luminosity) {
			yy[k] = TEFF->get(ODE->get_y(1,k))*(g/2.28e14) * 4.0*PI*1e10*R*R / (ZZ*ZZ);
			yy[k] = Lscale*yy[k] + (1.0-Lscale)*Lmin;
		} else {
			yy[k]=1.38e-16*pow((TEFF->get(ODE->get_y(1,k))*(g/2.28e14))/5.67e-5,0.25)/(1.6e-12*ZZ);
		}
	}
	TE.minit(xx,yy,nmodel);
	free_vector(xx,1,nmodel);
	free_vector(yy,1,nmodel);

	// calculate chisq
	double chisq=0.0;
	for (int i=1; i<=data.n; i++) {
		chisq += pow((data.TT[i] - TE.get(data.t[i]))/data.Te[i],2.0);	
		//printf("%lg %lg %lg\n", data.t[i], data.TT[i], TE.get(data.t[i]));
	}
	TE.tidy();
	printf("chisq = %lg\n", chisq);
	printf("chisq_nu = %lg/(%d-3) = %lg\n", chisq, data.n, chisq/(data.n-3));
	//return chisq;

}


