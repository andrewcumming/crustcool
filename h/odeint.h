class Crust;

class Ode_Int {
public:
  int ignore, kount, stiff, verbose, tri;
  double dxsav, minstep, hmax;
  void init(int n);
  void tidy(void);
  void go(double x1, double x2, double xstep,
	  double eps, Crust *crust);
  void set_bc(int n, double num);
  double get_x(int i);
  double get_y(int n, int i);
  double get_d(int n, int i);
  double *xp, **yp;
  int nok, nbad;

private:
  double **dydxp,*hstr,*ystart;
  int kmax,nvar;
  void rkck(double y[], double dydx[], int n, double x, double h,
	    double yout[],
	    double yerr[], Crust *crust);
  void rkqs(double y[], double dydx[], int n, double *x, double htry, 
	    double eps,	double yscal[], double *hdid, double *hnext,
	    Crust *crust);
  void odeint(double ystart[], int nvar, double x1, double x2, double eps, 
	      double h1,double hmin, int *nok, int *nbad,
	      Crust *crust);
#define float double
  void rk4(float y[], float dydx[], int n, float x, float h, float yout[],
	     Crust *crust);
  void rkdumb(float vstart[], int nvar, float x1, float x2, int nstep,
	Crust *crust);
  void rkscale(float vstart[], int nvar, float x1, float x2, float h1,
	Crust *crust);
#undef float 

#define float double
  float **d,*x;   // from stifbs.c
  
  void simpr(float y[], float dydx[], float dfdx[], float **dfdy, int n,
	     float xs, float htot, int nstep, float yout[],
	     Crust *crust);
  void bansimpr(float y[], float dydx[], float dfdx[], float **dfdy, int n,
	     float xs, float htot, int nstep, float yout[],
	     Crust *crust);
  void trisimpr(float y[], float dydx[], float dfdx[], float **dfdy, int n,
	     float xs, float htot, int nstep, float yout[],
	     Crust *crust);
  void stifbs(float y[], float dydx[], int nv, float *xx, float htry, float eps,
	      float yscal[], float *hdid, float *hnext,
		Crust *crust);
  void pzextr(int iest, float xest, float yest[], float yz[], float dy[], int nv);
  void lubksb(float **a, int n, int *indx, float b[]);
  void ludcmp(float **a, int n, int *indx, float *d);

  void tridag(float a[], float b[], float c[], float r[], float u[],
	      unsigned long n);
  void bandec(float **a, unsigned long n, int m1, int m2, float **al,
	      int *indx, float *d);
  void banbks(float **a, unsigned long n, int m1, int m2, float **al,
	      int *indx, float b[]);
#undef float
};

