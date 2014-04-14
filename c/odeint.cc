// odeint.cc
//
// 11th Mar 2002 -- added stiff integrator
//                  to use set "stiff" to 1
//

#include "../h/vector.h"
#include <stdio.h>
#include "math.h"
#include "../h/spline.h"
#include "../h/eos.h"
#include "../h/odeint.h"


void Ode_Int::tidy(void)
{
  this->nvar++;
  free_vector(this->ystart);
  free_vector(this->hstr);
  free_vector(this->xp);
  free_matrix(this->yp,this->nvar,this->kmax);
  free_matrix(this->dydxp,this->nvar,this->kmax);
}

void Ode_Int::init(int n, Ode_Int_Delegate *delegate)
{
	this->delegate = delegate;

  this->kmax=10000;
  this->nvar=n+1;
  this->ignore=0;
  this->dxsav=0.0;
  this->minstep=0.0;
	this->hmax=1e12;

  this->xp=vector(this->kmax);
  this->yp=matrix(this->nvar,this->kmax);
  this->dydxp=matrix(this->nvar,this->kmax);
  this->hstr=vector(this->kmax);
  this->ystart=vector(this->nvar);

  this->stiff = 0; // default is non-stiff eqns.
  this->verbose = 0; // turn off output
  this->tri=0; // don't assume a tridiagonal Jacobian

  this->nvar--;
}

void Ode_Int::set_bc(int n, double num)
{
  this->ystart[n]=num;
}

double Ode_Int::get_d(int n, int i)
{
  return this->dydxp[n][i];
}

double Ode_Int::get_x(int i)
{
  return this->xp[i];
}

double Ode_Int::get_y(int n, int i)
{
  return this->yp[n][i];
}

void Ode_Int::go(double x1, double x2, double xstep,
	double eps)
{
  if (this->dxsav == 0.0) this->dxsav=xstep;

  odeint(this->ystart,this->nvar,x1,x2,eps,xstep,this->minstep,&this->nok,
	 &this->nbad);
}

void Ode_Int::go_simple(double x1, double x2, int nstep)
{
  rkdumb(this->ystart,this->nvar,x1,x2,nstep);
  this->kount=nstep+1;
}


static double maxarg1,maxarg2;
#define FMAX(a,b) (maxarg1=(a),maxarg2=(b),(maxarg1) > (maxarg2) ?\
        (maxarg1) : (maxarg2))
static double minarg1,minarg2;
#define FMIN(a,b) (minarg1=(a),minarg2=(b),(minarg1) < (minarg2) ?\
        (minarg1) : (minarg2))
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))


#define SAFETY 0.9
#define PGROW -0.2
#define PSHRNK -0.25
#define ERRCON 1.89e-4
#define MAXSTP 10000000
#define TINY 1.0e-30

void Ode_Int::odeint(double ystart[], int nvar, double x1, double x2, 
		     double eps, double h1, double hmin, int *nok, int *nbad)
  // see NR pp 721.
{
  int nstp,i;
  double xsav=0.0,x,hnext,hdid,h;
  double *yscal,*y,*dydx;

  yscal=vector(nvar);
  y=vector(nvar);
  dydx=vector(nvar);
  x=x1;
  h=SIGN(h1,x2-x1);
  *nok = (*nbad) = kount =0;
  for (i=1;i<=nvar;i++) y[i]=ystart[i];

  if (kmax>0) xsav=x-this->dxsav*2.0;
  for (nstp=1;nstp<=MAXSTP;nstp++) {
	 if (this->verbose) printf("--- t=%lg delt=%lg nstep=%d\n", x, h,nstp);
    delegate->derivs(x,y,dydx);
    
    for (i=1;i<=nvar;i++)
      yscal[i]=fabs(y[i])+fabs(dydx[i]*h)+TINY;
      //yscal[i]=fabs(y[i])+TINY; 

    if (kmax > 0 && kount < kmax-1 && fabs(x-xsav) > fabs(this->dxsav)) {
      this->xp[++kount]=x; this->hstr[kount]=h;
      for(i=1;i<=nvar;i++) {
	this->yp[i][kount]=y[i];
	this->dydxp[i][kount]=dydx[i];
	//  printf("++++ %lg %lg %lg\n", this->xp[i], y[i], dydx[i]);
      }
      xsav=x;
    }
    if ((x+h-x2)*(x+h-x1) > 0.0) h=x2-x;
    //printf("h=%lg\n", h);
    if (this->stiff) {
      // if we've set the stiff flag use stifbs
      stifbs(y,dydx,nvar,&x,h,eps,yscal,&hdid,&hnext);
    } else {
      // otherwise rkqs
      rkqs(y,dydx,nvar,&x,h,eps,yscal,&hdid,&hnext);
    }
    //  printf("hdid=%lg, hnext=%lg\n", hdid, hnext);
    if (hdid == h) ++(*nok); else ++(*nbad);
    if ((x-x2)*(x2-x1) >= 0.0) {
      for (i=1;i<=nvar;i++) ystart[i]=y[i];
      if (kmax) {
	this->xp[++kount]=x;
	for (i=1;i<=nvar;i++) {
	  this->yp[i][kount]=y[i];
	  this->dydxp[i][kount]=dydx[i];
	  //	  printf("--++++ %lg %lg %lg\n", this->xp[i], y[i], dydx[i]);
	}      
      }
      free_vector(dydx);
      free_vector(y);
      free_vector(yscal);
      return;
    }
    if (fabs(hnext) <= hmin) printf("Step size too small in odeint..\n");
	//printf("%lg %lg %lg\n", this->hmax, h, hnext);
	if (hnext < this->hmax) h=hnext; else h=this->hmax;
  }
  printf("Too many steps in odeint...\n");
}


void Ode_Int::rkqs(double y[], double dydx[], int n, double *x, 
		   double htry, double eps, double yscal[], double *hdid, 
		   double *hnext)

  // Numerical recipes pp 719 

{
   int i;
   double errmax,h,htemp,xnew,*yerr,*ytemp;

   yerr=vector(n);
   ytemp=vector(n);
   h=htry;
   for (;;) {
      rkck(y,dydx,n,*x,h,ytemp,yerr);
      errmax=0.0;
      for (i=1;i<=(n-this->ignore);i++) {
	//printf("%lg %lg %lg %lg\n", y[i],  yerr[i], yscal[i], yerr[i]/yscal[i]);
	errmax=FMAX(errmax,fabs(yerr[i]/yscal[i]));
      }
      /*    for (i=1;i<=n;i++) errmax=FMAX(errmax,fabs(yerr[i]/yscal[i]));*/
      errmax /= eps;
      if (errmax <= 1.0) break;
      htemp=SAFETY*h*pow(errmax,PSHRNK);
      h=(h >= 0.0 ? FMAX(htemp,0.1*h) : FMIN(htemp,0.1*h));
      xnew=(*x)+h;
      if (xnew == *x) printf("stepsize underflow in rkqs!!!");
   }
   if (errmax > ERRCON) *hnext=SAFETY*h*pow(errmax,PGROW);
   else *hnext=5.0*h;
   *x += (*hdid=h);
   for (i=1;i<=n;i++) y[i]=ytemp[i];
   free_vector(ytemp);
   free_vector(yerr);
}


void Ode_Int::rkck(double y[], double dydx[], int n, double x, 
		   double h, double yout[], double yerr[])
{
   int i;
   static double a2=0.2,a3=0.3,a4=0.6,a5=1.0,a6=0.875,b21=0.2,
      b31=3.0/40.0,b32=9.0/40.0,b41=0.3,b42=-0.9,b43=1.2,
      b51=-11.0/54.0,b52=2.5,b53=-70.0/27.0,b54=35.0/27.0,
      b61=1631.0/55296.0,b62=175.0/512.0,b63=575.0/13824.0,
      b64=44275.0/110592.0,b65=253.0/4096.0,c1=37.0/378.0,
      c3=250.0/621.0,c4=125.0/594.0,c6=512.0/1771.0,
      dc5=-277.0/14336.0;
   double dc1=c1-2825.0/27648.0,dc3=c3-18575.0/48384.0,
      dc4=c4-13525.0/55296.0,dc6=c6-0.25;
   double *ak2,*ak3,*ak4,*ak5,*ak6,*ytemp;
   ak2=vector(n);
   ak3=vector(n);
   ak4=vector(n);
   ak5=vector(n);
   ak6=vector(n);
   ytemp=vector(n);
   for (i=1;i<=n;i++)
      ytemp[i]=y[i]+b21*h*dydx[i];
   delegate->derivs(x+a2*h,ytemp,ak2);
   for (i=1;i<=n;i++)
      ytemp[i]=y[i]+h*(b31*dydx[i]+b32*ak2[i]);
   delegate->derivs(x+a3*h,ytemp,ak3);
   for (i=1;i<=n;i++)
      ytemp[i]=y[i]+h*(b41*dydx[i]+b42*ak2[i]+b43*ak3[i]);
   delegate->derivs(x+a4*h,ytemp,ak4);
   for (i=1;i<=n;i++)
     ytemp[i]=y[i]+h*(b51*dydx[i]+b52*ak2[i]+b53*ak3[i]+b54*ak4[i]);
   delegate->derivs(x+a5*h,ytemp,ak5);
   for (i=1;i<=n;i++)
     ytemp[i]=y[i]+h*(b61*dydx[i]+b62*ak2[i]+b63*ak3[i]+b64*ak4[i]+b65*ak5[i]);
   delegate->derivs(x+a6*h,ytemp,ak6);
   for (i=1;i<=n;i++)
      yout[i]=y[i]+h*(c1*dydx[i]+c3*ak3[i]+c4*ak4[i]+c6*ak6[i]);
   for (i=1;i<=n;i++) {
     //printf("&&&&& %lg %lg %lg %lg %lg\n", yerr[i], ak3[i], ak4[i], ak5[i], ak6[i]);
      yerr[i]=h*(dc1*dydx[i]+dc3*ak3[i]+dc4*ak4[i]+dc5*ak5[i]+dc6*ak6[i]);
   }
   free_vector(ytemp);
   free_vector(ak6);
   free_vector(ak5);
   free_vector(ak4);
   free_vector(ak3);
   free_vector(ak2);
}



void Ode_Int::rkscale(double vstart[], int nvar, double x1, double x2, double h1)
{
	int i,k;
	double x,h;
	double *v,*vout,*dv;
	double hsum;

	v=vector(nvar);
	vout=vector(nvar);
	dv=vector(nvar);
	for (i=1;i<=nvar;i++) {
		v[i]=vstart[i];
		this->yp[i][1]=v[i];
	}
	this->xp[1]=x1;
	x=x1;
	h=h1; // user supplies initial step size
	k=1; // count the number of steps
	while ( x2>(x+h) ) {
		delegate->derivs(x,v,dv);

		// choose a step size
		hsum=0.0;
		for (i=1;i<=nvar;i++) {
		  hsum+=fabs(dv[i]/v[i]);
		}
		h=0.1/hsum;
		//printf("x=%lg, h=%lg\n  ", x, h);

		rk4(v,dv,nvar,x,h,vout);
		if ((double)(x+h) == x) printf("Step size too small in routine rkdumb\n");
		x += h;
		this->xp[k+1]=x;
		for (i=1;i<=nvar;i++) {
			v[i]=vout[i];
			this->yp[i][k+1]=v[i];
			this->dydxp[i][k+1]=dv[i];
		}

		// choose a new step size
		//hsum=0.0;
		//for (i=1;i<=nvar;i++) {
		//  hsum+=fabs(dv[i]/v[i]);
		//	}
		//h=0.1/hsum;
		//printf("x=%lg, h=%lg\n  ", x, h);
		k++;
	}
	this->kount=k;
	free_vector(dv);
	free_vector(vout);
	free_vector(v);
}

void Ode_Int::rkdumb(double vstart[], int nvar, double x1, double x2, int nstep)
{
	int i,k;
	double x,h;
	double *v,*vout,*dv;

	v=vector(nvar);
	vout=vector(nvar);
	dv=vector(nvar);
	for (i=1;i<=nvar;i++) {
		v[i]=vstart[i];
		this->yp[i][1]=v[i];
	}
	this->xp[1]=x1;
	x=x1;
	h=(x2-x1)/nstep;
	for (k=1;k<=nstep;k++) {
		delegate->derivs(x,v,dv);
		rk4(v,dv,nvar,x,h,vout);
		if ((double)(x+h) == x) printf("Step size too small in routine rkdumb\n");
		x += h;
		this->xp[k+1]=x;
		for (i=1;i<=nvar;i++) {
			v[i]=vout[i];
			this->yp[i][k+1]=v[i];
		}
	}
	free_vector(dv);
	free_vector(vout);
	free_vector(v);
}


void Ode_Int::rk4(double y[], double dydx[], int n, double x, double h, double yout[])
{
	int i;
	double xh,hh,h6,*dym,*dyt,*yt;

	dym=vector(n);
	dyt=vector(n);
	yt=vector(n);
	hh=h*0.5;
	h6=h/6.0;
	xh=x+hh;
	for (i=1;i<=n;i++) yt[i]=y[i]+hh*dydx[i];
	delegate->derivs(xh,yt,dyt);
	for (i=1;i<=n;i++) yt[i]=y[i]+hh*dyt[i];
	delegate->derivs(xh,yt,dym);
	for (i=1;i<=n;i++) {
		yt[i]=y[i]+h*dym[i];
		dym[i] += dyt[i];
	}
	delegate->derivs(x+h,yt,dyt);
	for (i=1;i<=n;i++)
		yout[i]=y[i]+h6*(dydx[i]+dyt[i]+2.0*dym[i]);
	free_vector(yt);
	free_vector(dyt);
	free_vector(dym);
}



// ------------------------ stiff integration routines -----------------------------



void Ode_Int::simpr(double y[], double dydx[], double dfdx[], double **dfdy, int n,
	double xs, double htot, int nstep, double yout[])
{
  //	void lubksb(double **a, int n, int *indx, double b[]);
  //	void ludcmp(double **a, int n, int *indx, double *d);
	int i,j,nn,*indx;
	double d,h,x,**a,*del,*ytemp;

	indx = new int [n+1];
	a=matrix(n,n);
	del=vector(n);
	ytemp=vector(n);
	h=htot/nstep;
	for (i=1;i<=n;i++) {
		for (j=1;j<=n;j++) a[i][j] = -h*dfdy[i][j];
		++a[i][i];
	}
	ludcmp(a,n,indx,&d);
	for (i=1;i<=n;i++)
		yout[i]=h*(dydx[i]+h*dfdx[i]);
	lubksb(a,n,indx,yout);
	for (i=1;i<=n;i++)
		ytemp[i]=y[i]+(del[i]=yout[i]);
	x=xs+h;
	delegate->derivs(x,ytemp,yout);
	for (nn=2;nn<=nstep;nn++) {
		for (i=1;i<=n;i++)
			yout[i]=h*yout[i]-del[i];
		lubksb(a,n,indx,yout);
		for (i=1;i<=n;i++)
			ytemp[i] += (del[i] += 2.0*yout[i]);
		x += h;
		delegate->derivs(x,ytemp,yout);
	}
	for (i=1;i<=n;i++)
		yout[i]=h*yout[i]-del[i];
	lubksb(a,n,indx,yout);
	for (i=1;i<=n;i++)
		yout[i] += ytemp[i];
	free_vector(ytemp);
	free_vector(del);
	free_matrix(a,n,n);
	delete [] indx;
	
}

void Ode_Int::bansimpr(double y[], double dydx[], double dfdx[], double **dfdy, int n,
	double xs, double htot, int nstep, double yout[])
{
	int i,j,nn, *indx;
	double d,h,x,**a,**al,*del,*ytemp;

	indx = new int [n+1];
	a=matrix(n,7);
	al=matrix(n,3);
	del=vector(n);
	ytemp=vector(n);
	h=htot/nstep;
	// fill up the compact representation of the band diagonal matrix a
	for (i=1; i<=n; i++) {
	  for (j=i-3; j<=i+3; j++) {
	    if (j>0) {
	      a[i][j-i+4] = -h*dfdy[i][j];
	    } else {
	      a[i][j-i+4] = 0.0;
	    }
	  }
	  a[i][4]+=1.0;   // diagonal element
	}
	bandec(a,n,3,3,al,indx,&d);    // replaces   ludcmp(a,n,indx,&d);
	for (i=1;i<=n;i++)
		yout[i]=h*(dydx[i]+h*dfdx[i]);
	banbks(a,n,3,3,al,indx,yout);    // here and later replaces   lubksb(a,n,indx,yout);
	for (i=1;i<=n;i++)
		ytemp[i]=y[i]+(del[i]=yout[i]);
	x=xs+h;
	delegate->derivs(x,ytemp,yout);
	for (nn=2;nn<=nstep;nn++) {
		for (i=1;i<=n;i++)
			yout[i]=h*yout[i]-del[i];
		banbks(a,n,3,3,al,indx,yout);
		for (i=1;i<=n;i++)
			ytemp[i] += (del[i] += 2.0*yout[i]);
		x += h;
		delegate->derivs(x,ytemp,yout);
	}
	for (i=1;i<=n;i++)
		yout[i]=h*yout[i]-del[i];
	banbks(a,n,3,3,al,indx,yout);
	for (i=1;i<=n;i++)
		yout[i] += ytemp[i];
	free_vector(ytemp);
	free_vector(del);
	free_matrix(a,n,n);
	delete [] indx;
}

void Ode_Int::trisimpr(double y[], double dydx[], double dfdx[], double **dfdy, int n,
	double xs, double htot, int nstep, double yout[])
{
	int i,nn;
	double h,x,*del,*ytemp;
	double *r,*AA,*BB,*CC;

	del=vector(n);
	ytemp=vector(n);
	r=vector(n);
	AA=vector(n);
	BB=vector(n);
	CC=vector(n);
	h=htot/nstep;
	for (i=1; i<=n; i++) BB[i]=1.0-h*dfdy[i][i];
	for (i=2; i<=n; i++) AA[i]=-h*dfdy[i][i-1];
	for (i=1; i<n; i++) CC[i]=-h*dfdy[i][i+1];

	for (i=1;i<=n;i++) r[i]=h*(dydx[i]+h*dfdx[i]);
	tridag(AA,BB,CC,r,yout,n);

	for (i=1;i<=n;i++)
		ytemp[i]=y[i]+(del[i]=yout[i]);
	x=xs+h;
	delegate->derivs(x,ytemp,yout);
	for (nn=2;nn<=nstep;nn++) {
		for (i=1;i<=n;i++)
		  r[i]=h*yout[i]-del[i];
		tridag(AA,BB,CC,r,yout,n);
		for (i=1;i<=n;i++)
			ytemp[i] += (del[i] += 2.0*yout[i]);
		x += h;
		delegate->derivs(x,ytemp,yout);
	}
	for (i=1;i<=n;i++) r[i]=h*yout[i]-del[i];
	tridag(AA,BB,CC,r,yout,n);
	for (i=1;i<=n;i++)
		yout[i] += ytemp[i];
	free_vector(ytemp);
	free_vector(del);

	free_vector(r);
	free_vector(AA);
	free_vector(BB);
	free_vector(CC);
}

#define KMAXX 7
#define IMAXX (KMAXX+1)
#define SAFE1 0.25
#define SAFE2 0.7
#define REDMAX 1.0e-5
#define REDMIN 0.7
#define TINY 1.0e-30
#define SCALMX 0.1


void Ode_Int::stifbs(double y[], double dydx[], int nv, double *xx, double htry, double eps,
	double yscal[], double *hdid, double *hnext)
{

	int i,iq,k,kk,km=0;
	static int first=1,kmax,kopt,nvold = -1;
	static double epsold = -1.0,xnew;
	double eps1,errmax,fact,h,red=1.0,scale=1.0,work,wrkmin,xest;
	double *dfdx,**dfdy,*err,*yerr,*ysav,*yseq;
	static double a[IMAXX+1];
	static double alf[KMAXX+1][KMAXX+1];
	static int nseq[IMAXX+1]={0,2,6,10,14,22,34,50,70};
	int reduct,exitflag=0;

	d=matrix(nv,KMAXX);
	dfdx=vector(nv);
	dfdy=matrix(nv,nv);
	err=vector(KMAXX);
	x=vector(KMAXX);
	yerr=vector(nv);
	ysav=vector(nv);
	yseq=vector(nv);
	for (int i=1; i<=nv; i++) {
	  dfdx[i]=0.0;
	  for (int j=1; j<=nv; j++)
	    dfdy[i][j]=0.0;
	}

	if(eps != epsold || nv != nvold) {
		*hnext = xnew = -1.0e29;
		eps1=SAFE1*eps;
		a[1]=nseq[1]+1;
		for (k=1;k<=KMAXX;k++) a[k+1]=a[k]+nseq[k+1];
		for (iq=2;iq<=KMAXX;iq++) {
			for (k=1;k<iq;k++)
				alf[k][iq]=pow(eps1,((a[k+1]-a[iq+1])/
					((a[iq+1]-a[1]+1.0)*(2*k+1))));
		}
		epsold=eps;
		nvold=nv;
		a[1] += nv;
		for (k=1;k<=KMAXX;k++) a[k+1]=a[k]+nseq[k+1];
		for (kopt=2;kopt<KMAXX;kopt++)
			if (a[kopt+1] > a[kopt]*alf[kopt-1][kopt]) break;
		kmax=kopt;
	}
	h=htry;
	for (i=1;i<=nv;i++) ysav[i]=y[i];
	delegate->jacobn(*xx,y,dfdx,dfdy,nv);
	if (*xx != xnew || h != (*hnext)) {
		first=1;
		kopt=kmax;
	}
	reduct=0;
	for (;;) {
		for (k=1;k<=kmax;k++) {
			xnew=(*xx)+h;
			if (xnew == (*xx)) { printf("stepsize underflow in stifbs, h=%lg\n",h); exitflag=1;} 
			if (this->tri==1) trisimpr(ysav,dydx,dfdx,dfdy,nv,*xx,h,nseq[k],yseq);
			else if (this->tri==2) bansimpr(ysav,dydx,dfdx,dfdy,nv,*xx,h,nseq[k],yseq);
			else simpr(ysav,dydx,dfdx,dfdy,nv,*xx,h,nseq[k],yseq);
			double sqrarg = (h/nseq[k]);
			xest=(sqrarg == 0.0 ? 0.0 : sqrarg*sqrarg);
			pzextr(k,xest,yseq,y,yerr,nv);
			if (k != 1) {
				errmax=TINY;
				for (i=1;i<=nv;i++) errmax=FMAX(errmax,fabs(yerr[i]/yscal[i]));
				errmax /= eps;
				km=k-1;
				err[km]=pow(errmax/SAFE1,1.0/(2*km+1));
			}
			if (k != 1 && (k >= kopt-1 || first)) {
				if (errmax < 1.0) {
					exitflag=1;
					break;
				}
				if (k == kmax || k == kopt+1) {
					red=SAFE2/err[km];
					break;
				}
				else if (k == kopt && alf[kopt-1][kopt] < err[km]) {
						red=1.0/err[km];
						break;
					}
				else if (kopt == kmax && alf[km][kmax-1] < err[km]) {
						red=alf[km][kmax-1]*SAFE2/err[km];
						break;
					}
				else if (alf[km][kopt] < err[km]) {
					red=alf[km][kopt-1]/err[km];
					break;
				}
			}
		}
		if (exitflag) break;
		red=FMIN(red,REDMIN);
		red=FMAX(red,REDMAX);
		h *= red;
		reduct=1;
	}
	*xx=xnew;
	*hdid=h;
	first=0;
	wrkmin=1.0e35;
	for (kk=1;kk<=km;kk++) {
		fact=FMAX(err[kk],SCALMX);
		work=fact*a[kk+1];
		if (work < wrkmin) {
			scale=fact;
			wrkmin=work;
			kopt=kk+1;
		}
	}
	*hnext=h/scale;
	if (kopt >= k && kopt != kmax && !reduct) {
		fact=FMAX(scale/alf[kopt-1][kopt],SCALMX);
		if (a[kopt+1]*fact <= wrkmin) {
			*hnext=h/fact;
			kopt++;
		}
	}
	free_vector(yseq);
	free_vector(ysav);
	free_vector(yerr);
	free_vector(x);
	free_vector(err);
	free_matrix(dfdy,nv,nv);
	free_vector(dfdx);
	free_matrix(d,nv,KMAXX);
}

#undef KMAXX
#undef IMAXX
#undef SAFE1
#undef SAFE2
#undef REDMAX
#undef REDMIN
#undef TINY
#undef SCALMX





void Ode_Int::pzextr(int iest, double xest, double yest[], double yz[], double dy[], int nv)
{
	int k1,j;
	double q,f2,f1,delta,*c;

	c=vector(nv);
	x[iest]=xest;
	for (j=1;j<=nv;j++) dy[j]=yz[j]=yest[j];
	if (iest == 1) {
		for (j=1;j<=nv;j++) d[j][1]=yest[j];
	} else {
		for (j=1;j<=nv;j++) c[j]=yest[j];
		for (k1=1;k1<iest;k1++) {
			delta=1.0/(x[iest-k1]-xest);
			f1=xest*delta;
			f2=x[iest-k1]*delta;
			for (j=1;j<=nv;j++) {
				q=d[j][k1];
				d[j][k1]=dy[j];
				delta=c[j]-q;
				dy[j]=f1*delta;
				c[j]=f2*delta;
				yz[j] += dy[j];
			}
		}
		for (j=1;j<=nv;j++) d[j][iest]=dy[j];
	}
	free_vector(c);
}





void Ode_Int::lubksb(double **a, int n, int *indx, double b[])
{
	int i,ii=0,ip,j;
	double sum;

	for (i=1;i<=n;i++) {
		ip=indx[i];
		sum=b[ip];
		b[ip]=b[i];
		if (ii)
			for (j=ii;j<=i-1;j++) sum -= a[i][j]*b[j];
		else if (sum) ii=i;
		b[i]=sum;
	}
	for (i=n;i>=1;i--) {
		sum=b[i];
		for (j=i+1;j<=n;j++) sum -= a[i][j]*b[j];
		b[i]=sum/a[i][i];
	}
}



#define TINY 1.0e-20;

void Ode_Int::ludcmp(double **a, int n, int *indx, double *d)
{
	int i,imax=0,j,k;
	double big,dum,sum,temp;
	double *vv;

	vv=vector(n);
	*d=1.0;
	for (i=1;i<=n;i++) {
		big=0.0;
		for (j=1;j<=n;j++)
			if ((temp=fabs(a[i][j])) > big) big=temp;
		if (big == 0.0) printf("Singular matrix in routine ludcmp\n");
		vv[i]=1.0/big;
	}
	for (j=1;j<=n;j++) {
		for (i=1;i<j;i++) {
			sum=a[i][j];
			for (k=1;k<i;k++) sum -= a[i][k]*a[k][j];
			a[i][j]=sum;
		}
		big=0.0;
		for (i=j;i<=n;i++) {
			sum=a[i][j];
			for (k=1;k<j;k++)
				sum -= a[i][k]*a[k][j];
			a[i][j]=sum;
			if ( (dum=vv[i]*fabs(sum)) >= big) {
				big=dum;
				imax=i;
			}
		}
		if (j != imax) {
			for (k=1;k<=n;k++) {
				dum=a[imax][k];
				a[imax][k]=a[j][k];
				a[j][k]=dum;
			}
			*d = -(*d);
			vv[imax]=vv[j];
		}
		indx[j]=imax;
		if (a[j][j] == 0.0) a[j][j]=TINY;
		if (j != n) {
			dum=1.0/(a[j][j]);
			for (i=j+1;i<=n;i++) a[i][j] *= dum;
		}
	}
	free_vector(vv);
}
#undef TINY



/* note #undef's at end of file */


void Ode_Int::tridag(double a[], double b[], double c[], double r[], double u[],
	unsigned long n)
{
	unsigned long j;
	double bet,*gam;

	gam=vector(n);
	if (b[1] == 0.0) printf("Error 1 in tridag\n");
	u[1]=r[1]/(bet=b[1]);
	for (j=2;j<=n;j++) {
		gam[j]=c[j-1]/bet;
		bet=b[j]-a[j]*gam[j];
		if (bet == 0.0)	printf("Error 2 in tridag\n");
		u[j]=(r[j]-a[j]*u[j-1])/bet;
	}
	for (j=(n-1);j>=1;j--)
		u[j] -= gam[j+1]*u[j+1];
	free_vector(gam);
}



#define SWAP(a,b) {dum=(a);(a)=(b);(b)=dum;}
#define TINY 1.0e-20

void Ode_Int::bandec(double **a, unsigned long n, int m1, int m2, double **al,
	int *indx, double *d)
{
	unsigned long i,j,k,l;
	int mm;
	double dum;

	mm=m1+m2+1;
	l=m1;
	for (i=1;i<=m1;i++) {
		for (j=m1+2-i;j<=mm;j++) a[i][j-l]=a[i][j];
		l--;
		for (j=mm-l;j<=mm;j++) a[i][j]=0.0;
	}
	*d=1.0;
	l=m1;
	for (k=1;k<=n;k++) {
		dum=a[k][1];
		i=k;
		if (l < n) l++;
		for (j=k+1;j<=l;j++) {
			if (fabs(a[j][1]) > fabs(dum)) {
				dum=a[j][1];
				i=j;
			}
		}
		indx[k]=i;
		if (dum == 0.0) a[k][1]=TINY;
		if (i != k) {
			*d = -(*d);
			for (j=1;j<=mm;j++) SWAP(a[k][j],a[i][j])
		}
		for (i=k+1;i<=l;i++) {
			dum=a[i][1]/a[k][1];
			al[k][i-k]=dum;
			for (j=2;j<=mm;j++) a[i][j-1]=a[i][j]-dum*a[k][j];
			a[i][mm]=0.0;
		}
	}
}
#undef SWAP
#undef TINY


/* note #undef's at end of file */
#define SWAP(a,b) {dum=(a);(a)=(b);(b)=dum;}

void Ode_Int::banbks(double **a, unsigned long n, int m1, int m2, double **al,
	int *indx, double b[])
{
	unsigned long i,k,l;
	int mm;
	double dum;

	mm=m1+m2+1;
	l=m1;
	for (k=1;k<=n;k++) {
		i=indx[k];
		if (i != k) SWAP(b[k],b[i])
		if (l < n) l++;
		for (i=k+1;i<=l;i++) b[i] -= al[k][i-k]*b[k];
	}
	l=1;
	for (i=n;i>=1;i--) {
		dum=b[i];
		for (k=2;k<=l;k++) dum -= a[i][k]*b[k+i-1];
		b[i]=dum/a[i][1];
		if (l < mm) l++;
	}
}
#undef SWAP


