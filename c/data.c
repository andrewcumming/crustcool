#include <stdio.h>
#include <math.h>
#include "../h/spline.h"
#include "../h/nr.h"
#include "../h/nrutil.h"
#include "../h/odeint.h"

struct Data {
	double *t, *TT, *Te;
	int n;
} data;
	
void read_in_data(const char *fname) 
{
	if (1) {   



			/*
			// hardcode the data for 1731
			double t0=51930.5;
			data.n=7;    // 7 data points for BC09
			data.t = vector(1,data.n);
			data.TT = vector(1,data.n);
			data.Te = vector(1,data.n);

			data.t[1]=51995.1; data.TT[1]=103.2; data.Te[1]=1.7;
			data.t[2]=52165.7; data.TT[2]=88.9; data.Te[2]=1.3;
			data.t[3]=52681.6; data.TT[3]=75.5; data.Te[3]=2.2;
			data.t[4]=52859.5; data.TT[4]=73.3; data.Te[4]=2.3;
			data.t[5]=53430.5; data.TT[5]=71.0; data.Te[5]=1.8;
			data.t[6]=53500.4; data.TT[6]=66.0; data.Te[6]=4.5;
			data.t[7]=53525.4; data.TT[7]=70.3; data.Te[7]=2.1;
			data.t[8]=54969.7; data.TT[8]=63.1; data.Te[8]=2.1;

			for (int i=1; i<=data.n; i++) {
				data.t[i]-=t0;
			}
			*/

			/*
			// hardcode the data for 1659
			double t0=52159.5;
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

			for (int i=1; i<=data.n; i++) {
				data.t[i]-=t0;
			}
			*/


		// hardcode the data for XTEJ
		double t0=0.0;
		data.n=13;
		data.t = vector(1,data.n);
		data.TT = vector(1,data.n);
		data.Te = vector(1,data.n);
		data.t[1]=3.12; data.TT[1]=165.387; data.Te[1]=3.776;
		data.t[2]=10.98; data.TT[2]=161.689; data.Te[2]=2.710;
		data.t[3]=16.39; data.TT[3]=156.146; data.Te[3]=1.392;
		data.t[4]=49.66; data.TT[4]=150.495; data.Te[4]=1.218;
		data.t[5]=174.50; data.TT[5]=138.791; data.Te[5]=4.009;
		data.t[6]=298.47; data.TT[6]=138.366; data.Te[6]=1.892;
		data.t[7]=431.24; data.TT[7]=131.826; data.Te[7]=2.817;
		data.t[8]=540.25; data.TT[8]=127.959; data.Te[8]=1.578;
		data.t[9]=592.85; data.TT[9]=134.298; data.Te[9]=2.251;
		data.t[10]=652.79; data.TT[10]=128.455; data.Te[10]=2.100;
		data.t[11]=705.55; data.TT[11]=127.682; data.Te[11]=1.848;
		data.t[12]=795.80; data.TT[12]=127.785; data.Te[12]=1.657;
		data.t[13]=1159.01; data.TT[13]=124.719; data.Te[13]=1.740;
		//	data.t[14]=1905.96; data.TT[14]=112.500; data.Te[14]=3.423;
		// remove last data point to compare to Page & Reddy
	
		for (int i=1; i<=data.n; i++) data.t[i]-=t0;
	
	} else {	
	
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



double calculate_chisq(Ode_Int *ODE, Spline *TEFF, double g, double ZZ)	
// uses the result of the cooling to calculate chi-squared
{
	int nmodel = ODE->kount;
	
	// set up a spline which has the prediction for observed Teff vs time 
	Spline TE;
	double *yy = vector(1,nmodel);
	double *xx = vector(1,nmodel);
	for (int k=1; k<=nmodel; k++) { 
		xx[k]=ODE->get_x(k)*ZZ/(3600.0*24.0);
		yy[k]=1.38e-16*pow((TEFF->get(ODE->get_y(1,k))*(g/2.28e14))/5.67e-5,0.25)/(1.6e-12*ZZ);
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
	return chisq;

}


