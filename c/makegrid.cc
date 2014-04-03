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
	
	printf("Enter log10 base column of He layer  (0 to force iron)..."); scanf("%lg",&yi);
	printf("Enter B field in G (0 for unmagnetized)..."); scanf("%lg",&Bfield);

	Envelope envelope;
	envelope.use_potek_eos_in_He=0;
	envelope.use_potek_cond_in_He=0;
	envelope.use_potek_eos_in_Fe=0;
	envelope.use_potek_cond_in_Fe=0;
	if (Bfield > 0.0) envelope.use_potek_kff=1;
	else envelope.use_potek_kff=0;
	envelope.make_grid(yi,Bfield);   // results are in "out/grid"
}
