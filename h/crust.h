class Ode_Int;
class Eos;
class Spline;

class Crust {
public:
	Crust();
       ~Crust();
	void setup(void);
	void evolve(double time, double mdot);
	
	int N;
	double Pb, Pt, yt, dx;
	double *P, *rho, *T;

	int output, use_my_envelope, gpe;
	
	double mass,radius,g,ZZ;
	
	double timesofar, last_time_output;
	double Lmin, Lscale, Qimp, B, kncrit;
	int gap,accr,use_potek_eos;		
			
	int force_precalc,extra_heating,nuflag,accreting,force_cooling_bc;
	double rhot,rhob,heating_P1,heating_P2;
	double energy_deposited_outer,energy_deposited_inner,energy_slope;
	double mdot,outburst_duration;
	double extra_y,extra_Q,deep_heating_factor;
	
	double angle_mu,Tt,Tc,Qrho,Qinner;
	
	Spline TEFF;
	
	Ode_Int ODE;
	
	void derivs(double t, double T[], double dTdt[]);
	void jacobn(double, double *, double *, double **, int);
					
private:
	int hardwireQ;
	double *CP, *K, *F, *NU, *EPS, *Qheat, *Qimpur;

	int nbeta;
	double **CP_grid, **K1_grid, **K0_grid, **NU_grid, **EPS_grid, **KAPPA_grid, **K1perp_grid, **K0perp_grid;
	double betamin, betamax, deltabeta;
	FILE *fp,*fp2;
	
	void set_up_grid(const char *fname);
	void get_TbTeff_relation(void);
	void set_composition(void);
	double crust_heating(int i);
	double energy_deposited(int i);
	void precalculate_vars(void);
	
	void output_result_for_step(int j, FILE *fp, FILE *fp2,double timesofar,double *last_time_output);
	
	double dTdt(int i, double *T);
	void calculate_vars(int i, double T, double y, double *CP, double *K, double *NU, double *EPS);
	void outer_boundary(double T1, double K1, double CP1, double NU1, double EPS1, double *T0, double *K0, double *CP0, double *NU0, double *EPS0);
	void inner_boundary(double TN, double KN, double CPN, double NUN, double EPSN, double *TN1, double *KN1, double *CPN1, double *NUN1, double *EPSN1);
	double calculate_heat_flux(int i, double *T);
	
	Spline AASpline; 
	Spline ZZSpline;
	Spline YnSpline;
	Eos *EOS;
};
