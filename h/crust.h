#include "../h/odeint.h"
#include "../h/spline.h"
#include "../h/eos.h"

struct GridPoint {
	double rho, CP, P, K, F, T, NU, EPS, Qheat, Qimpur, r;
};


class Crust: public Ode_Int_Delegate {
public:
	Crust();
    ~Crust();
	void setup(void);
	void evolve(double time, double mdot);
	void set_temperature_profile(double *rhovec,double *Tvec,int nvec);
	
	int N;
	double Pb, Pt, yt, dx;

	GridPoint *grid;

	int output, use_my_envelope, gpe, resume;
	
	double mass,radius,g,ZZ;
	double C_core, Lnu_core_norm, Lnu_core_alpha;
	
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

	int nbeta;
	double **CP_grid, **K1_grid, **K0_grid, **NU_grid, **EPS_grid, **KAPPA_grid, **K1perp_grid, **K0perp_grid;
	double betamin, betamax, deltabeta;
	FILE *fp,*fp2;
	
	void set_up_grid(const char *fname);
	void get_TbTeff_relation(void);
	void set_composition(void);
	double crust_heating(int i);
	void precalculate_vars(void);
	double eps_from_heat_source(double P,double y1,double y2,double Q_heat);

	void output_result_for_step(int j, FILE *fp, FILE *fp2,double timesofar,double *last_time_output);
	void read_T_profile_from_file(void);
	
	double dTdt(int i, double *T);
	void calculate_vars(int i);
	double calculate_heat_flux(int i, double *T);
	void outer_boundary(void);
	
	Spline AASpline; 
	Spline ZZSpline;
	Spline YnSpline;
	Eos *EOS;
};
