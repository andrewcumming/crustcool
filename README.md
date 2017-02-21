### Overview and setup

This code follows the thermal evolution of a neutron star crust. It is designed to model observations of accreting neutron stars in quiescence and the decline of magnetar outbursts.

To compile, `make` should compile the cooling code `crustcool`.

#### Requirements

* GNU Scientific Library [GSL](http://www.gnu.org/software/gsl/) 
* [condegin13.f](http://www.ioffe.ru/astro/conduct/condin.html) fortran routine by A. Potekhin to calculate thermal conductivity (put this in the `c` directory)
* [eos14.f](http://www.ioffe.ru/astro/EIP/index.html) and [eosmag14.f](http://www.ioffe.ru/astro/EIP/index.html) fortran routines by A. Potekhin which calculate equation of state (put this in the `c` directory)
* For MCMC, [emcee](https://github.com/dfm/emcee) and [triangle_plot](https://github.com/dfm/triangle.py)

### Setting parameters 

The file `init.dat` sets up the run. The parameters are

	Tt		temperature at the top of the crust during accretion
	Tc		core temperature
	Qimp	impurity parameter
	Qinner 	(optional) a different impurity parameter for the inner crust
	Qrho	the density at which the Q changes from Qimp to Qinner (default 1e12)

	cooling_bc	if set to 1, use a cooling b.c. at the top even during accretion,
				otherwise keep the temperature at the top fixed (to the value Tt) 
				during accretion.
	extra_heating	if set, turn on extra heating in the NS ocean during accretion.
	extra_Q		strength of extra heating in MeV
	extra_y		depth of extra heating in g/cm^2
	deep_heating_factor	factor by which to multiply the deep heating strength (useful to be able to control
				deep heating and shallow heating independently; equivalent to changing mdot for deep heating only)

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
	resume	start new output files with an intially isothermal crust (=0), or resume using the temperature profile from last time and append to the output(=1) 
	
	SFgap	neutron superfluid gap. Choices are
			0=normal neutrons (not SF)
			1=SFB03 (default choice), 2=AWPII, 3=Gaussian Tc(k)
			4=all or nothing, the neutrons are normal for k<kncrit
			5=B1 from Reddy and Page
			6=BCS from Reddy and Page
	kncrit	neutrons are normal for kn<kncrit (to use this set SFgap=4)
	
	potek_eos	use the EOS routines from Potekhin in the EOS

	piecewise	if =1 then the initial temperature is specified in a piecewise
				format in the lines beginning with > in this file
	timetorun	time to run in days
	neutrinos	include neutrino cooling (1=yes 0=no)

	toutburst	accretion outburst duration in years
	accreted	crust composition  1=accreted crust (HZ1990)
					0=equilibrium crust (HP1994;DH2001)
					2=accreted crust (HZ2003)
					
	C_core		core heat capacity at 1e8 K
	Lnu_core_norm	normalization in expression for core neutrino luminosity at 1e8 K
	Lnu_core_alpha	temperature sensitivity of core neutrino luminosity (e.g. slow process=8, fast process=6)
	
	
You can include comments (`#`) in the `init.dat` file, blank lines are ignored, and lines beginning with `>` are (optionally) to specify the piecewise initial temperature profile. A pair of double comment symbols `##` can be used to comment out a block of lines.

If you give an argument, e.g.

	crustcool source

then the code will look for the file `init/init.dat.source` instead of `init.dat`. This is useful to keep different setups for modelling different data sets for example.

### Example

	crustcool 1659_example
	python plot_tc.py

should give the following plot in `tc.pdf`

![MXB 1659 lightcurve](https://github.com/andrewcumming/crustcool/raw/master/1659_example.png)


### MCMC

`mcee.py` is an MCMC driver which uses the [emcee](http://dan.iel.fm/emcee/current) python code. 

`mcplot.py` plots the output of `mcee.py`. It uses the [triangle_plot](http://pypi.python.org/pypi/triangle_plot) plotting routines.

`mcmc.py` is a python driver for MCMC using a simple Metropolis algorithm.


### Published cooling curves from this code

* [Horowitz et al. (2015)](http://arxiv.org/abs/1410.2197) MXB 1659-29 (Figs 2 and 3)
* [Scholz et al. (2014)](http://lanl.arxiv.org/abs/1401.6965) Swift J1822.3-1606 (Fig. 2)
* [Cackett et al. (2013)](http://arxiv.org/abs/1306.1776) MXB 1659-29 (Figure 4)
* [An et al. (2013)](http://arxiv.org/abs/1212.0184) CXOU J164710.2-455216 (Fig. 1, see `init.dat.1647`)
* [Scholz et al. (2012)](http://lanl.arxiv.org/abs/1204.1034) Swift J1822.3-1606 (Fig. 8, see `init.dat.1822`)
* [An et al. (2012)](http://arxiv.org/abs/1208.1419) SGR 1627-41 (Fig.3, 1998 and 2008 outbursts, see `init.dat.1627`)
