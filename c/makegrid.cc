// makegrid.cc
//
// Makes envelope models
//
 
#include <stdio.h>
#include "../h/envelope.h"

//------------------------------------------------------------------------


int main(void)
{
	double yi, Bfield;
	
	printf("Enter log10 base column of He layer ..."); scanf("%lg",&yi);
	printf("Enter B field in G (0 for unmagnetized)..."); scanf("%lg",&Bfield);

	Envelope envelope;
	envelope.use_potek_eos_in_He=0;
	envelope.use_potek_eos_in_Fe=0;
	envelope.use_potek_cond_in_Fe=1;
	envelope.use_potek_cond_in_He=1;
	envelope.make_grid(yi,Bfield);   // results are in "out/grid"
}
