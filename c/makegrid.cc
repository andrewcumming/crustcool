// makegrid.cc
//
// The idea here is to calculate a range of constant flux solutions 
// on a constant grid in y, so that we get a mapping between 
// T_b,y_b and the flux
 
#include <stdio.h>
#include <string.h>
#include "math.h"
#include <stdarg.h>

#include "../h/nr.h"
#include "../h/nrutil.h"
#include "../h/envelope.h"

//------------------------------------------------------------------------


int main(void)
{
	double yi, Bfield;
	
	printf("Enter log10 base column of He layer ..."); scanf("%lg",&yi);
	printf("Enter B field in G (0 for unmagnetized)..."); scanf("%lg",&Bfield);

	Envelope envelope;
	envelope.make_grid(yi,Bfield);
}
