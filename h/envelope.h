#include "../h/odeint.h"
#include "../h/eos.h"
#include <stdio.h>

class Envelope: public Ode_Int_Delegate {
public:
	Envelope();
	~Envelope();
	
	double g;
	double F;
	double Bfield;
	double yi;
	FILE *fp;
	
	void derivs(double t, double T[], double dTdt[]);
    void jacobn(double, double *, double *, double **, int);

	void doint(double F);
	void calculate(void);

	void make_grid(double yi, double B);

private:
	static double Wrapper_find_surf_eqn(double r);
	double find_surf_eqn(double r);
	Eos *EOS;
	Ode_Int ODE, ODE2;
};
