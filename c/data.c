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

		if (!strncmp(sourcename,"1822",4)) {  //  Swift J1822

			data.luminosity=1;

			t0=55756.533180;
			data.n=59;
			data.t = vector(1,data.n);
			data.TT = vector(1,data.n);
			data.Te = vector(1,data.n);
			data.t[1]=55769.2; data.TT[1]=1.18534e-10; data.Te[1]=1.2384e-12;
			data.t[2]=55769.5; data.TT[2]=1.02749e-10; data.Te[2]=8.7312e-12;
			data.t[3]=55770.3; data.TT[3]=9.84182e-11; data.Te[3]=9.75336e-12;
			data.t[4]=55771.2; data.TT[4]=9.27598e-11; data.Te[4]=8.2469e-12;
			data.t[5]=55772.3; data.TT[5]=8.96419e-11; data.Te[5]=7.37617e-12;
			data.t[6]=55777.1; data.TT[6]=7.11526e-11; data.Te[6]=9.8682e-13;
			data.t[7]=55778; data.TT[7]=7.07113e-11; data.Te[7]=7.15515e-12;
			data.t[8]=55779; data.TT[8]=6.5864e-11; data.Te[8]=6.35643e-12;
			data.t[9]=55780.4; data.TT[9]=5.86284e-11; data.Te[9]=5.14376e-12;
			data.t[10]=55781.4; data.TT[10]=5.89166e-11; data.Te[10]=6.78554e-12;
			data.t[11]=55786.4; data.TT[11]=4.93738e-11; data.Te[11]=5.54141e-12;
			data.t[12]=55787.6; data.TT[12]=4.67744e-11; data.Te[12]=3.86155e-12;
			data.t[13]=55788.1; data.TT[13]=4.49513e-11; data.Te[13]=4.11844e-12;
			data.t[14]=55789.5; data.TT[14]=4.25396e-11; data.Te[14]=4.06365e-12;
			data.t[15]=55790.3; data.TT[15]=3.69914e-11; data.Te[15]=4.03279e-12;
			data.t[16]=55800.8; data.TT[16]=2.94367e-11; data.Te[16]=3.40334e-12;
			data.t[17]=55807.2; data.TT[17]=2.52772e-11; data.Te[17]=4.92286e-12;
			data.t[18]=55822.7; data.TT[18]=1.78974e-11; data.Te[18]=1.73994e-12;
			data.t[19]=55822.7; data.TT[19]=2.074e-11; data.Te[19]=6.7678e-13;
			data.t[20]=55824.5; data.TT[20]=2.02612e-11; data.Te[20]=3.84792e-12;
			data.t[21]=55829.1; data.TT[21]=2.06206e-11; data.Te[21]=3.00107e-12;
			data.t[22]=55835.1; data.TT[22]=1.75026e-11; data.Te[22]=3.17129e-12;
			data.t[23]=55841.7; data.TT[23]=1.48954e-11; data.Te[23]=2.40557e-12;
			data.t[24]=55849.2; data.TT[24]=1.27647e-11; data.Te[24]=2.53515e-12;
			data.t[25]=55856.4; data.TT[25]=1.17095e-11; data.Te[25]=3.29858e-12;
			data.t[26]=55862.1; data.TT[26]=9.10405e-12; data.Te[26]=9.6258e-13;
			data.t[27]=55867.1; data.TT[27]=9.90377e-12; data.Te[27]=2.96883e-13;
			data.t[28]=55976.4; data.TT[28]=5.17721e-12; data.Te[28]=9.41038e-13;
			data.t[29]=55977; data.TT[29]=5.4422e-12; data.Te[29]=8.86704e-13;
			data.t[30]=55978.1; data.TT[30]=4.99704e-12; data.Te[30]=9.60791e-13;
			data.t[31]=55981.9; data.TT[31]=4.98466e-12; data.Te[31]=9.5639e-13;
			data.t[32]=55982.9; data.TT[32]=4.79108e-12; data.Te[32]=9.93451e-13;
			data.t[33]=55985; data.TT[33]=4.94233e-12; data.Te[33]=8.54705e-13;
			data.t[34]=55991.1; data.TT[34]=4.84548e-12; data.Te[34]=9.21239e-13;
			data.t[35]=56031.1; data.TT[35]=4.62992e-12; data.Te[35]=1.00324e-12;
			data.t[36]=56037; data.TT[36]=4.53741e-12; data.Te[36]=1.69035e-13;
			data.t[37]=56052.6; data.TT[37]=4.48224e-12; data.Te[37]=8.58434e-13;
			data.t[38]=56073; data.TT[38]=4.2973e-12; data.Te[38]=9.30267e-13;
			data.t[39]=56095.5; data.TT[39]=3.76885e-12; data.Te[39]=9.20909e-13;
			data.t[40]=56104.1; data.TT[40]=4.56548e-12; data.Te[40]=9.71978e-13;
			data.t[41]=56114.2; data.TT[41]=3.78601e-12; data.Te[41]=9.09639e-13;
			data.t[42]=56156.1; data.TT[42]=3.4532e-12; data.Te[42]=8.81743e-13;
			data.t[43]=56161.5; data.TT[43]=3.42478e-12; data.Te[43]=9.46656e-13;
			data.t[44]=56206; data.TT[44]=3.00719e-12; data.Te[44]=1.20071e-12;
			data.t[45]=56238.5; data.TT[45]=2.28136e-12; data.Te[45]=8.4258e-13;
			data.t[46]=56334.7; data.TT[46]=1.82984e-12; data.Te[46]=5.62029e-13;
			data.t[47]=56355.1; data.TT[47]=1.84399e-12; data.Te[47]=5.61143e-13;
			data.t[48]=56386.8; data.TT[48]=1.85315e-12; data.Te[48]=1.00599e-12;
			data.t[49]=56387; data.TT[49]=1.92864e-12; data.Te[49]=6.3453e-13;
			data.t[50]=56408.6; data.TT[50]=1.43908e-12; data.Te[50]=6.22505e-13;
			data.t[51]=56456.1; data.TT[51]=1.12244e-12; data.Te[51]=4.664e-13;
			data.t[52]=56459.3; data.TT[52]=1.05311e-12; data.Te[52]=5.083e-13;
			data.t[53]=56490.4; data.TT[53]=1.02683e-12; data.Te[53]=7.54752e-13;
			data.t[54]=56491.1; data.TT[54]=1.04981e-12; data.Te[54]=3.03876e-13;
			data.t[55]=56535.3; data.TT[55]=7.95675e-13; data.Te[55]=3.79581e-13;
			data.t[56]=56598; data.TT[56]=5.11623e-13; data.Te[56]=6.49733e-13;
			data.t[57]=56599.1; data.TT[57]=4.95097e-13; data.Te[57]=3.90062e-13;
			data.t[58]=56713; data.TT[58]=1.316e-13; data.Te[58]=2.84e-14;
			data.t[59]=56714.8; data.TT[59]=1.4031e-13; data.Te[59]=2.87e-14;
			
			// convert fluxes to luminosity
			double distance=1.6;   //distance in kpc
			for (int i=1; i<=data.n; i++) {
				data.TT[i]*=4.0*PI*pow(3.086e21*distance,2.0);
				data.Te[i]*=4.0*PI*pow(3.086e21*distance,2.0);
			}
			
		}

		if (!strncmp(sourcename,"1731",4)) {  // KS 1731
			t0=51930.5;
			data.n=8;    // 7 data points for BC09
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

			//	data.t[9]=t0+5230.0; data.TT[9]=54.0; data.Te[9]=2.1;
			//	data.t[9]=t0+5230.0; data.TT[9]=63.1; data.Te[9]=2.1;

		}

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
		
		if (!strncmp(sourcename,"xtej",4)) {  // XTE J1701	
			t0=0.0;
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
		}
		
		if (!strncmp(sourcename,"0556",4)) {  
			t0=0.0;
			data.n=9;
			data.t = vector(1,data.n);
			data.TT = vector(1,data.n);
			data.Te = vector(1,data.n);
			data.t[1]=5.5; data.TT[1]=309.0; data.Te[1]=7.0;
			data.t[2]=16.1; data.TT[2]=308.0; data.Te[2]=3.0;
			data.t[3]=23.3; data.TT[3]=296.0; data.Te[3]=3.0;
			//data.t[4]=31.9; data.TT[4]=361.0; data.Te[4]=10.0;
			data.t[4]=51.1; data.TT[4]=276.0; data.Te[4]=2.0;
			//data.t[6]=85.4; data.TT[6]=317.0; data.Te[6]=5.0;
			data.t[5]=104.5; data.TT[5]=250.3; data.Te[5]=0.9;
			data.t[6]=134.9; data.TT[6]=241.4; data.Te[6]=0.9;
			data.t[7]=150.9; data.TT[7]=239.0; data.Te[7]=2.0;
			data.t[8]=291.8; data.TT[8]=207.0; data.Te[8]=2.0;
			data.t[9]=497.1; data.TT[9]=182.9; data.Te[9]=1.0;
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


