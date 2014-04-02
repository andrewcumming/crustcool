#include "../h/envelope.h"
#include "../h/nr.h"
#include "math.h"

void* pt2EnvObject;

Envelope::Envelope()
{	
	pt2EnvObject=(void*) this;
  	
	static Eos myeos(1);
    this->EOS = &myeos;
	this->EOS->X[1]=1.0;
	this->EOS->Qimp=0.0;
	this->EOS->accr=1;

	this->ODE.init(1);
	this->ODE2.init(1);
	
	// Default values
	this->g = 2.28e14;
	
}


Envelope::~Envelope() 
{
	this->ODE.tidy();
	this->ODE2.tidy();
}


void Envelope::make_grid(double yi, double B)
// yi is the log10 of the base column of the He layer
// B is the magnetic field strength
{
	this->fp = fopen("out/grid","w");
	this->yi = yi;
	this->Bfield = B;

	for (int i=0; i<=80; i++) {
		double F=19.0+i*0.1;
	   	this->doint(F);
		printf("."); fflush(stdout);
    	for (int j=1; j<=ODE2.kount; j++) {
			fprintf(this->fp, "%lg %lg %lg\n", this->ODE2.get_x(j), log10(this->ODE2.get_y(1,j)), F);
		}
    	for (int j=1; j<=ODE.kount; j++) {
			fprintf(this->fp, "%lg %lg %lg\n", this->ODE.get_x(j), log10(this->ODE.get_y(1,j)), F);
		}
	}
	printf("\n");
	fclose(this->fp);
}


void Envelope::doint(double F)
// for the specified flux, integrate inwards to see if we match the base temperature
{
  	F=pow(10.0,F);  // flux is provided as log_10
  	this->F=F;

  	// we do this in two steps: light element layer first
	this->EOS->B=this->Bfield;
	this->EOS->use_potek_eos = 0;
	this->EOS->use_potek_cond = 1;
	this->EOS->A[1]=4.0; this->EOS->Z[1]=2.0;   // Helium
	//this->EOS->A[1]=56.0; this->EOS->Z[1]=26.0;  // Iron

  	// set surface temperature. We integrate from tau=2/3
	//printf("----\n%g %g\n----\n",this->Wrapper_find_surf_eqn(1e-3),this->Wrapper_find_surf_eqn(1e5));
	double yt=1.0; //zbrent(this->Wrapper_find_surf_eqn,1e-3,1e5,1e-5);
  	if (yt==1e-5 || yt==1e3) printf("yt out of bounds (%lg)\n", yt);

  	double Tt=pow(F/5.67e-5,0.25);
	this->ODE2.set_bc(1,Tt);

  	// integrate
  	this->ODE2.go_simple(log10(yt),this->yi,(int)((this->yi-log10(yt))/0.05),dynamic_cast<Ode_Int_Delegate *>(this));

  	// keep the base temperature for the next integration
  	double base_T=this->ODE2.get_y(1,this->ODE2.kount);

  	// setup ocean
  	this->ODE.set_bc(1,base_T);
	this->EOS->Qimp=0.0;
	this->EOS->use_potek_cond = 1;
	this->EOS->use_potek_eos = 0;
 	this->EOS->A[1]=56.0; this->EOS->Z[1]=26.0;  // Iron
  
  	// integrate through the ocean to the desired depth 
	this->ODE.go_simple(this->yi,18.5,(int)((18.5-this->yi)/0.05),dynamic_cast<Ode_Int_Delegate *>(this));
}



double Envelope::Wrapper_find_surf_eqn(double y)
{
  Envelope* mySelf = (Envelope*) pt2EnvObject;
  return mySelf->find_surf_eqn(y);
}

double Envelope::find_surf_eqn(double y)
{
	this->EOS->P=this->g*y; this->EOS->T8=1e-8*pow(this->F/5.67e-5,0.25);
	this->EOS->rho=this->EOS->find_rho();
	//printf("%lg %lg %lg\n",this->EOS->P,this->EOS->rho,this->EOS->T8);
	return this->EOS->opac()*y-2.0/3.0;
}


void Envelope::derivs(double x, double ff[], double dfdx[])
// Evaluate derivative dT/dy for the integrator
{
  	double y=pow(10.0,x);
  	double T=ff[1];

  	// find density; the composition is already set
  	this->EOS->P=this->g*y; EOS->T8=T*1e-8; 
  	this->EOS->rho=EOS->find_rho();
	this->EOS->Yn=0.0;
	if (x > this->yi) EOS->set_composition_by_density();

  	double kappa=EOS->opac();
	if (0) {
		double mykcond = this->EOS->kcond;
		// here I call Potekhin's routine for the conductivity  
 		EOS->kcond=3.024e20*pow(EOS->T8,3)/(EOS->potek_cond()*EOS->rho);
		// and then swap it in
  		kappa=1.0/((1.0/EOS->kcond)+(1.0/EOS->kappa_rad));
	}

  	// Heat equation for constant flux
  	dfdx[1]=2.303*y*3305.1*this->F*kappa/pow(T,3.0);

	// if the temperature gradient is superadiabatic, set it to the adiabatic gradient (convection)
	double conv_grad = 2.303*T*EOS->del_ad();
	if (dfdx[1] > conv_grad) dfdx[1]=conv_grad;
}

void Envelope::jacobn(double, double *, double *, double **, int) 
{
}

