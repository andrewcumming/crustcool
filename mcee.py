#!/usr/bin/env python
import random
import numpy
import matplotlib
# matplotlib.use('macosx')
import matplotlib.pyplot as plt
import emcee
import uuid
import time
import subprocess
import re
import os
import datetime
import imessage
from math import sqrt, pi
import ns
import random
from multiprocessing import Pool

os.environ["OMP_NUM_THREADS"] = "1"
def main():

	nwalkers, ndim = 10, 3
	nsteps = 10
	# 200, 1000
	dir = '1659'
	if os.path.exists('mcmc/'+dir):
		print(('The output directory mcmc/'+dir+' already exists!'))
		exit()

	# clean up temporary files if there are any left over from a previous run
	os.system('rm /tmp/init.dat.*')
		
	# parameters are   x = [Tc7, Qimp, Tb8, mdot, M, R, kncrit]
	#p0 = emcee.utils.sample_ball([3.0,13.0],[0.5,0.5],nwalkers)
	#p0 = emcee.utils.sample_ball([3.0,13.0,1.4,10.0],[0.5,0.5,0.1,1.0],nwalkers)
	#p0 = emcee.utils.sample_ball([9.0,1.0,4.0,0.1,1.2,13.0],[0.3,0.1,0.3,0.01,0.01,0.5],nwalkers)
	#p0 = emcee.utils.sample_ball([7.7,0.5,4.0,0.1],[0.3,0.3,0.3,0.01],nwalkers)
	p0 = emcee.utils.sample_ball([7.7,0.5,4.0],[0.3,0.3,0.3],nwalkers)
	
	# Q,Lscale,Touter,Tc,M,R,rhotransition,Tinner
	#p0 = emcee.utils.sample_ball([1.0,0.2,9.0,2.0,1.8,2.0,10.0,1.0],
			#[0.1,0.03,0.1,0.1,0.07,0.3,0.1,0.2],nwalkers)
	# Q,Lscale,Edep,Tc,M,g
	#p0 = emcee.utils.sample_ball([1.0,0.2,3.0,2.0,1.8,2.0],[0.1,0.03,0.1,0.1,0.07,0.3],nwalkers)
	# Q,Lscale,Edep,Tc,M,R	
	#p0 = emcee.utils.sample_ball([1.0,0.2,3.0,2.0,1.4,15.0],[0.1,0.03,0.1,0.1,0.07,0.1],nwalkers)
	#p0 = emcee.utils.sample_ball([1.0,0.2,3.0,2.0,2.0,10.0],[0.1,0.03,0.1,0.1,0.07,0.3],nwalkers)
	#p0 = emcee.utils.sample_ball([0.0,0.2,3.0,2.0,1.6,12.0],[0.3,0.03,0.1,0.1,0.1,0.8],nwalkers)

	# for i in range(0,nwalkers):
	# 	p0[i][4] = 1.1 + random.random()*(2.4-1.1)
	# 	(g1,g2) = gravity_range(p0[i][4])
	# 	p0[i][5]=g1 + random.random()*(g2-g1)

	with Pool(8) as pool:
		sampler=emcee.EnsembleSampler(nwalkers,ndim,lnprob,pool=pool)

		time_per_step = (2411.0/4000.0)*1.2/6.0   # time to run one cooling curve
		print(('Starting run at ', str(datetime.datetime.now())))
		start_time = time.time()
		print(('Estimated time to run is ', time_per_step*nsteps*nwalkers, ' seconds'))
		print('Estimated completion time is ', str(datetime.datetime.now()+datetime.timedelta(seconds=time_per_step*nsteps*nwalkers)))
		sampler.run_mcmc(p0, nsteps, progress=True)
		print('time to run = ',time.time() - start_time,'seconds')
		print(("Mean acceptance fraction: {0:.3f}".format(numpy.mean(sampler.acceptance_fraction))))

# throw away the first burn steps as burn-in
	nburn = 0
	samples=sampler.chain[:,nburn:,:].reshape((-1,ndim))

	os.makedirs('mcmc/'+dir)
	outputFile = open('mcmc/'+dir+'/samples.dat','wb')
	numpy.save(outputFile, numpy.array(samples))
	outputFile.close()
	outputFile = open('mcmc/'+dir+'/sampler_chain.dat','wb')
	numpy.save(outputFile, numpy.array(sampler.chain))
	outputFile.close()
	outputFile = open('mcmc/'+dir+'/sampler_lnprobability.dat','wb')
	numpy.save(outputFile, numpy.array(sampler.lnprobability))
	outputFile.close()

	fp = open('mcmc/'+dir+'/chi.dat','w')
	probs=sampler.lnprobability.reshape([-1])
	for i,prob in enumerate(probs):
		print(i,-2.0*prob," ".join(str(x) for x in samples[i,:]), file=fp)
	fp.close()

	os.system('cp mcee.py mcmc/'+dir+'/')

	#imessage.send('mcee '+dir+' has finished running')


def set_params(x,name):
	data="""resume 0

mass	1.62
radius	11.2
timetorun	10000.0
toutburst 2.5
mdot	0.1

precalc 0
ngrid	50
SFgap	1
kncrit	0
accreted 1

Qimp	3.2
Tc	3.1e7

# fixed top temperature during accretion:
cooling_bc	0
extra_heating	0
Tt	4.2e8

# or cooling boundary condition with shallow heating:
#cooling_bc	1
#extra_heating	1
#extra_Q	1.4
#extra_y	1e13
"""
	# Q,Lscale,Edep,Tc,M,R	
	#radius = ns.R(x[4],x[5])
	#Lmin = PYL(x[3]*1e7, 1e14, x[4], radius)
	
	params = {
		'Tc':x[0]*1e7,
		'Qimp':10.0**x[1],
#		'Lscale':x[1],
#		'Edep':x[2],
		'Tt':x[2]*1e8,
#		'mdot':x[3],
#		'mass':x[4],
#		'radius':radius,
#		'Lmin':Lmin
#		'extra_Q':x[0],
#		'extra_y':x[1]
	}

	for key,value in list(params.items()):
		data = re.sub("%s(\\t?\\w?).*\\n" % key,"%s\\t%g\\n" % (key,value),data)

	fout = open('/tmp/init.dat.'+name,'w')
	fout.write(data)
	
	# for piecewise, x[2]=top temperature, x[6]=log10 transition density, 
	#x[7]=inner crust temperature
#	print('>0\t%g' % (x[2]*1e8), file=fout)
#	print('>%g\t%g\t%g' % (10.0**x[6],x[2]*1e8,x[7]*1e8), file=fout)
#	print('>1e14\t%g\t-1' % (x[7]*1e8), file=fout)
#	print('>-1\t-1', file=fout)
	
	fout.close()


def get_chisq(x):
	name = str(uuid.uuid4())
	set_params(x,name)
	# give crustcool a second parameter so that it looks in /tmp for the init.dat file
	data = subprocess.check_output( ["./crustcool",name,'1'])
	# print(data)
	chisq = float(re.search(b'chisq = ([-+]?[0-9]*\.?[0-9]+)',data).group(1))
#	print x[0],x[1],x[2],x[3],x[4],x[5],x[6],chisq
	os.system('rm /tmp/init.dat.'+name)
	return chisq


def lnprob(x):
	# minimum and maxiumum allowed values
	# (assume a flat prior within this range)

	# Q,Lscale,Edep,Tc,M,R
	#(g1,g2) = gravity_range(x[4])
#	xmin=numpy.array([-2.0,0.01,0.1,1.0,1.1,g1])
#	xmax=numpy.array([2.0,0.5,100.0,3.0,2.4,g2])
	# Q,Lscale,Touter,Tc,M,R,rhotransition,Tinner
#	xmin=numpy.array([-2.0,0.01,5.0,1.0,1.1,g1,9.0,0.1])
#	xmax=numpy.array([2.0,0.5,20.0,3.0,2.4,g2,12.0,3.0])
#	xmin=numpy.array([1.0,-3.0,0.1,0.0,1.1,8.0])
#	xmax=numpy.array([100.0,3.0,100.0,3.0,2.5,16.0])
	xmin=numpy.array([0.0,-3.0,0.0])
	xmax=numpy.array([100.0,3.0,100.0])
	#xmin=numpy.array([0.0,-3.0,0.0,0.0])
	#xmax=numpy.array([100.0,3.0,100.0,1.0])
	if (len((x<xmin).nonzero()[0])>0 or len((x>xmax).nonzero()[0])>0):
		return -numpy.inf
	# causal limit  radius has to be R > 4.36km (M/Msun)
	#if (ns.R(x[4],x[5]) < 4.36*x[4]):
	#	return -numpy.inf
	#print 'Trying ',x
	chisq=get_chisq(x)
	#print x[0],x[1],x[2],chisq
	return -chisq/2.0


def gravity_range(mass):
	# returns the allowed gravity range at a given mass
	g1 = ns.grav(mass,16.0)
	if g1<0.6:
		g1=0.6
	g2 = ns.grav(mass,8.0)
	if g2>5.0:   # limit gravity to be <5 (this keeps up causal even at 2.4 M_sun)
		g2=5.0
	return (g1,g2)

def PYL(Tc, B, mass, radius):
		# Calculates the expected luminosity from a dipole field using 
		# the Potekhin & Yakovlev envelope relations
		ZZ = 1.0/(1.0 - 2.0*6.67e-8*2e33*mass/(radius*1e5*9e20))**0.5
		grav = 6.67e-8*2e33*mass*ZZ/(radius**2*1e10)
		T9 = Tc/1e9
		xi = T9 - 0.001*(1e-14*grav)**0.25*sqrt(7.0*T9)
		flux = 5.67e-5 * 1e24 * grav*1e-14 * ((7*xi)**2.25+(0.333*xi)**1.25)
		# now correct for B ... 
		# use the enhancement along the field direction;
		# chi = 1.0 + 0.0492 * (1d-12*B)^0.292 / T9^0.24
		# flux *= chi^4
		# or use eq. (31) or PY2001  which gives F(B)/F(0)
		bet = 0.074*sqrt(1e-12*B)*T9**(-0.45)
		a1=5059.0*T9**0.75/sqrt(1.0 + 20.4*sqrt(T9) + 138.0*T9**1.5 + 1102.0*T9*T9);
		a2=1484.0*T9**0.75/sqrt(1.0 + 90.0*T9**1.5+ 125.0*T9*T9);
		a3=5530.0*T9**0.75/sqrt(1.0 + 8.16*sqrt(T9) + 107.8*T9**1.5+ 560.0*T9*T9);
		fac = (1.0 + a1*bet**2 + a2*bet**3 + 0.007*a3*bet**4)/(1.0+a3*bet**2);
		flux *= fac
		flux *= 4.0*pi*1e10*radius**2
		flux/=ZZ**2
		return flux


if __name__ == '__main__':
	main()
