#include "../h/odeint.h"
#include "../h/eos.h"
#include <stdio.h>

class Envelope: public Ode_Int_Delegate {
public:
	Envelope();
	~Envelope();
	
	double g;    // gravity in cgs units
	
	int use_potek_eos_in_He;
	int use_potek_eos_in_Fe;
	int use_potek_cond_in_He;
	int use_potek_cond_in_Fe;
	
	// makes a grid of envelope models ('out/grid')
	// yi is the He layer column depth (log10)
	// set B=0 for unmagnetized envelope	
	void make_grid(double yi, double B);

private:
	static double Wrapper_find_surf_eqn(double r);
	double find_surf_eqn(double r);
	Eos *EOS;
	Ode_Int ODE, ODE2;
	void doint(void);
	void calculate(void);
	FILE *fp;
	double Bfield;
	double yi;
	double F;	
	void derivs(double t, double T[], double dTdt[]);
    void jacobn(double, double *, double *, double **, int);
};
