#include <stdio.h>
#include <math.h>
#include <string.h>
#include "../h/crust.h"
#include "../h/data.h"

void Data::read_in_data(const char *sourcename) 
{
	printf("\nComparing with data: source = %s\n",sourcename);
	this->luminosity = 0;
	
	if (1) {    // hardcoded data    

		double t0=0.0;

		if (!strncmp(sourcename,"1659",4)) {  // MXB 1659	
			t0=52159.5;
			this->n=8;    // 7 data points in BC09
			this->t = new double[this->n+1];
			this->TT = new double[this->n+1];
			this->Te = new double[this->n+1];

			this->t[1]=52197.8; this->TT[1]= 121; this->Te[1]= 1;
			this->t[2]=52563.2; this->TT[2]= 85; this->Te[2]= 1;
			this->t[3]=52712.2; this->TT[3]= 77; this->Te[3]= 1;
			this->t[4]=52768.9; this->TT[4]= 73; this->Te[4]= 1;
			this->t[5]=53560.0; this->TT[5]= 58; this->Te[5]= 2;
			this->t[6]=53576.7; this->TT[6]= 54; this->Te[6]= 3;
			this->t[7]=54583.8; this->TT[7]= 56; this->Te[7]= 2;
			this->t[8]=56113; this->TT[8]= 48.8; this->Te[8]= 1.6;
		}
		
		for (int i=1; i<=this->n; i++) {
			this->t[i]-=t0;
		}
	
	} else {	
	
		char fname[100]="data/";
		strcat(fname,sourcename);
	
		printf("Reading data from %s\n", fname);
	
		FILE *fp = fopen(fname,"r");
	
		double t0;
		fscanf(fp, "%lg %d\n", &t0, &this->n);
	
		this->t = new double[this->n+1];
		this->TT = new double[this->n+1];
		this->Te = new double[this->n+1];
		
		for (int i=1; i<=this->n; i++) {
			double dummy;
			fscanf(fp,"%lg %lg %lg %lg %lg\n",  &this->t[i], &dummy, &dummy, &this->TT[i],&this->Te[i]);
			this->t[i]-=t0;
		}	
		fclose(fp);	
	}
}



void Data::calculate_chisq(Crust &crust)	
//void Data::calculate_chisq(Ode_Int *ODE, Spline *TEFF, double g, double ZZ, double R,double Lscale,double Lmin)	
// uses the result of the cooling to calculate chi-squared
{
	int nmodel = crust.ODE.kount;
	
	double g=crust.g;
	double ZZ=crust.ZZ;
	double R=crust.radius;
	double Lscale=crust.Lscale;
	double Lmin=crust.Lmin;
	
	// set up a spline which has the prediction for observed Teff vs time 
	Spline TE;
	double *yy = new double[nmodel+1];
	double *xx = new double[nmodel+1];
	for (int k=1; k<=nmodel; k++) { 
		xx[k]=crust.ODE.get_x(k)*ZZ/(3600.0*24.0);
		if (this->luminosity) {
			yy[k] = crust.TEFF.get(crust.ODE.get_y(1,k))*(g/2.28e14) * 4.0*M_PI*1e10*R*R / (ZZ*ZZ);
			yy[k] = Lscale*yy[k] + (1.0-Lscale)*Lmin;
		} else {
			yy[k]=1.38e-16*pow((crust.TEFF.get(crust.ODE.get_y(1,k))*(g/2.28e14))/5.67e-5,0.25)/(1.6e-12*ZZ);
		}
	}
	TE.minit(xx,yy,nmodel);
	delete [] xx;
	delete [] yy;

	// calculate chisq
	double chisq=0.0;
	for (int i=1; i<=this->n; i++) {
		chisq += pow((this->TT[i] - TE.get(this->t[i]))/this->Te[i],2.0);	
		//printf("%lg %lg %lg\n", this->t[i], this->TT[i], TE.get(this->t[i]));
	}
	TE.tidy();
	printf("chisq = %lg\n", chisq);
	printf("chisq_nu = %lg/(%d-3) = %lg\n", chisq, this->n, chisq/(this->n-3));
}


