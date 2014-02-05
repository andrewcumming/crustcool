## Overview and setup

This code follows the thermal evolution of a neutron star crust. It is designed to model observations of accreting neutron stars in quiescence and the decline of magnetar outbursts.

To compile, `make` should compile the cooling code `crustcool`. You will need to create the directory `o` for the object files, and a directory `gon_out` for output during the runs.

## Setting parameters 

The file `init.dat` sets up the run. The parameters are

	Tt		temperature at the top of the crust during accretion
	Tc		core temperature
	Qimp	impurity parameter
	Qinner 	(optional) a different impurity parameter for the inner crust
	Qrho	the density at which the Q changes from Qimp to Qinner (default 1e12)

	latent_heat	include (=1) or don't include (=0) the latent heat
	convection	include (=1) or don't include (=0) a simple model of convective fluxes f
						from compositionally-driven convection in the liquid layer (Medin & Cumming 2012)
	cooling_bc	if set to 1, use a cooling b.c. at the top even during accretion,
				otherwise keep the temperature at the top fixed (to the value Tt) 
				during accretion.
	extra_heating	if set, turn on extra heating in the NS ocean during accretion.
	extra_Q		strength of extra heating in MeV
	extra_y		depth of extra heating in g/cm^2

	Edep	energy deposited (used for magnetar heating)
	Einner	(optional) a different value of energy deposited for the inner crust
	rhot	lowest density to be heated
	rhob	highest density to be heated
	energy_slope	energy deposited is multiplied by  (rho_10)**energy_slope

	mass	neutron star mass in solar masses
	radius	neutron star radius in km

	gpe		1=iron envelope (use out/grid_He4 as the outer boundary)
			0=He envelope (use out/grid_He9 as the outer boundary; same as BC09)

	Bfield  magnetic field strength in the crust in G
	angle_mu	determines the Teff-Tb relation used for B>0. If angle_mu = -1 (default)
		then an angle-averaged relation is used; otherwise angle_mu in the range 0 to 1
		specifies the local angle of the field relative to the vertical

	envelope	if =1 then use the calculated envelop (out/grid_..) for B>0; if =0, then use
				the analytic envelope from the literature

	mdot	accretion rate in Eddington units (1.0 == 8.8e4 g/cm^2/s)

	precalc	force a precalc (1) or instead load in previously saved precalc (0)
	ngrid	number of grid points
	ytop	column depth at the top of the grid (default 1e12)
	output	write output files (=1) or suppress output (=0) (e.g. for mcmc we don't need output)
	
	SFgap	neutron superfluid gap. Choices are
			0=normal neutrons (not SF)
			1=SFB03, 2=AWPII, 3=Gaussian Tc(k)
			4=all or nothing, the neutrons are normal for k<kncrit
			5=B1 from Reddy and Page
			6=BCS from Reddy and Page
	kncrit	neutrons are normal for kn<kncrit (to use this set SFgap=4)
	sph		whether to inlcude SF phonons (0=no 1=yes)

	potek_eos	use the EOS routines from Potekhin in the EOS

	piecewise	if =1 then the initial temperature is specified in a piecewise
				format in the lines beginning with > in this file
	timetorun	time to run in days
	neutrinos	include neutrino cooling (1=yes 0=no)
	instant		heat "instantly" if =1, otherwise model the outburst

	toutburst	accretion outburst duration in years
	accreted	crust composition  1=accreted crust (HZ1990)
					0=equilibrium crust (HP1994;DH2001)
					2=accreted crust (HZ2003)
	
You can include comments (`#`) in the `init.dat` file, blank lines are ignored, and lines beginning with `>` are (optionally) to specify the piecewise initial temperature profile. A pair of double comment symbols `##` can be used to comment out a block of lines.

If you give an argument, e.g.

	crustcool source

then the code will look for the file `init/init.dat.source` instead of `init.dat`. This is useful to keep different setups for modelling different data sets for example.

### MCMC

`mcmc.py` is a python driver for MCMC using a simple Metropolis algorithm. For each set of trial parameters, it creates an `init.dat` file based on `init/init.dat.mcmc`. The output of `crustcool` is redirected to a file `tmp` which is then scanned for the chi-squared value.

`mcee.py` is an MCMC driver which uses the [emcee](http://dan.iel.fm/emcee/current) python code. 

`mcplot.py` plots the output of `mcee.py`. It uses the [triangle_plot](http://pypi.python.org/pypi/triangle_plot) plotting routines.


### Published cooling curves from this code

* An et al. (2013) Fig. 1 (see `init.dat.1647`)
* Scholz et al. (2012) Fig. 8 (`init.dat.1822`)
* Scholz et al. 2014 submitted paper on 1822
* An et al. (2012) Fig.3 (1998 and 2008 outbursts, see `init.dat.1627`)
