class Ode_Int_Delegate {
public:
	virtual void derivs(double t, double T[], double dTdt[]){};
	virtual void jacobn(double, double *, double *, double **, int){};
};


class Ode_Int {
public:
  int ignore, kount, stiff, verbose, tri;
  double dxsav, minstep, hmax;
  void init(int n,Ode_Int_Delegate *delegate);
  void tidy(void);
  void go(double x1, double x2, double xstep, double eps);
  void go_simple(double x1, double x2, int nstep);	
  void set_bc(int n, double num);
  double get_x(int i);
  double get_y(int n, int i);
  double get_d(int n, int i);
  double *xp, **yp;
  int nok, nbad;
  Ode_Int_Delegate *delegate;

private:
  double **dydxp,*hstr,*ystart;
  int kmax,nvar;
  void rkck(double y[], double dydx[], int n, double x, double h,
	    double yout[],
	    double yerr[]);
  void rkqs(double y[], double dydx[], int n, double *x, double htry, 
	    double eps,	double yscal[], double *hdid, double *hnext);
  void odeint(double ystart[], int nvar, double x1, double x2, double eps, 
	      double h1,double hmin, int *nok, int *nbad);
  void rk4(double y[], double dydx[], int n, double x, double h, double yout[]);
  void rkdumb(double vstart[], int nvar, double x1, double x2, int nstep);
	void rkscale(double vstart[], int nvar, double x1, double x2, double h1);

  double **d,*x;   // from stifbs.c
  
  void simpr(double y[], double dydx[], double dfdx[], double **dfdy, int n,
	     double xs, double htot, int nstep, double yout[]);
  void bansimpr(double y[], double dydx[], double dfdx[], double **dfdy, int n,
	     double xs, double htot, int nstep, double yout[]);
  void trisimpr(double y[], double dydx[], double dfdx[], double **dfdy, int n,
	     double xs, double htot, int nstep, double yout[]);
  void stifbs(double y[], double dydx[], int nv, double *xx, double htry, double eps,
	      double yscal[], double *hdid, double *hnext);
  void pzextr(int iest, double xest, double yest[], double yz[], double dy[], int nv);
  void lubksb(double **a, int n, int *indx, double b[]);
  void ludcmp(double **a, int n, int *indx, double *d);

  void tridag(double a[], double b[], double c[], double r[], double u[],
	      unsigned long n);
  void bandec(double **a, unsigned long n, int m1, int m2, double **al,
	      int *indx, double *d);
  void banbks(double **a, unsigned long n, int m1, int m2, double **al,
	      int *indx, double b[]);
};


static double maxarg1,maxarg2;
#define FMAX(a,b) (maxarg1=(a),maxarg2=(b),(maxarg1) > (maxarg2) ?\
        (maxarg1) : (maxarg2))
static double minarg1,minarg2;
#define FMIN(a,b) (minarg1=(a),minarg2=(b),(minarg1) < (minarg2) ?\
        (minarg1) : (minarg2))
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))



