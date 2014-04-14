// These functions to declare vectors and matrices
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>

double *vector(long nl, long nh)
{
	double *v = new double [nh+1];
	return v;
}

double **matrix(long nrl, long nrh, long ncl, long nch)
{
	double **m = new double*[nrh+1];
	for (int i = 0; i < nrh+1; ++i)
		m[i] = new double[nch+1];		
	return m;
}

void free_vector(double *v, long nl, long nh)
{
	delete [] v;
}

void free_matrix(double **m, long nrl, long nrh, long ncl, long nch)
{
	for (int i = 0; i<nrh+1; ++i) delete [] m[i];
	delete [] m;
}



