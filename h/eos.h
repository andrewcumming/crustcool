class Eos {
public:
	Eos(int n);
	~Eos();
	// density, temperature and composition
	double P, rho, T8, *A, *Z, *X, Qimp, Yn;

	// other parameters
	double gamma_melt, kncrit, B;
	int accr, gap;

	// find density from pressure
	double find_rho(void);
	static double Wrapper_find_rho_eqn(double r);
	double find_rho_eqn(double r);

	// mean molecular weights
	double Ye(void);
	double Yi(void);
	double YZ2(void);
	void set_composition_by_pressure(void);
	void set_composition_by_density(void);
	double set_Ye, set_Yi, set_YZ2;

	// equation of state
	double pe(void);
	double pemod(void);
	double ptot(void);
	double Utot(void);
	double Chabrier_EF(void);
	double Fermi_Inv_1_2(double F);
	double FermiI(int k, double T8, double EF);
	double x(void);
	double gamma(void);
	double Uex(void);
	double eta(void);

	// thermodynamics
	double f(void);
	double CP(void);
	double CV(void);
	double cvion, cv_alpha, cvrad, cve, cvneut;
	double del_ad(void);
	double chi(double *x);
	double Gamma1(void);

	// conductivity and neutrinos
	double eps_nu(void);
	double TC(void);
	double Q1,Q2,Q3,Q4,Q5,Q6;
	double K_cond(double ef);
	double f_ee, f_ei, f_eQ, f_ep;
	double econd(void);
	double opac(void);
	double kes, kff, kcond, kgaunt, kappa_rad;
	double gff(double Z, double eta);
	double J(double x,double y);
	double lambda2;
	double lamei(int n);
	double Fep(int flag);
	double Eep(double q);

	// interface to Potekhin's routines
	int use_potek_cond, use_potek_eos;
	double potek_cond(void);
	void potek_eos(double *P_out, double *cv_out_i, double *cv_out_e);
	double Kperp;

private:
	int ns;  // number of species
	double Fermi_n, Fermi_alpha;
	double expint(int n, double x);	

};
