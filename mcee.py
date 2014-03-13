import random
import numpy
import matplotlib
matplotlib.use('macosx')
import matplotlib.pyplot as plt
import emcee
import uuid
import time
import subprocess
import re
import os
import datetime
from math import sqrt, pi


def main():

	nwalkers, ndim = 200, 6
	dir = '1822_5'
	if os.path.exists('mcmc/'+dir):
		print 'The output directory mcmc/'+dir+' already exists!'
		exit()

	# clean up temporary files if there are any left over from a previous run
	os.system('rm /tmp/init.dat.*')
		
	# parameters are   x = [Tc7, Qimp, Tb8, mdot, M, R, kncrit]
	#p0 = emcee.utils.sample_ball([3.0,13.0],[0.5,0.5],nwalkers)
	#p0 = emcee.utils.sample_ball([3.0,13.0,1.4,10.0],[0.5,0.5,0.1,1.0],nwalkers)
	#p0 = emcee.utils.sample_ball([9.0,1.0,4.0,0.1,1.2,13.0],[0.3,0.1,0.3,0.01,0.01,0.5],nwalkers)
	#p0 = emcee.utils.sample_ball([7.7,0.5,4.0,0.1],[0.3,0.3,0.3,0.01],nwalkers)
	#p0 = emcee.utils.sample_ball([7.7,0.5,4.0],[0.3,0.3,0.3],nwalkers)
	
	# Q,Lscale,Edep,Tc,M,R
	p0 = emcee.utils.sample_ball([0.0,0.2,3.0,2.0,1.6,12.0],[0.3,0.03,0.1,0.1,0.1,0.8],nwalkers)

	sampler=emcee.EnsembleSampler(nwalkers,ndim,lnprob,threads=20)

	print 'Starting run at ', str(datetime.datetime.now())
	start_time = time.time()
	sampler.run_mcmc(p0, 500)
	print 'time to run = ',time.time() - start_time,'seconds'
	print("Mean acceptance fraction: {0:.3f}"
	                .format(numpy.mean(sampler.acceptance_fraction)))

	# throw away the first burn steps as burn-in
	nburn = 0
	samples=sampler.chain[:,nburn:,:].reshape((-1,ndim))

	os.makedirs('mcmc/'+dir)
	outputFile = open('mcmc/'+dir+'/samples.dat','w')
	numpy.save(outputFile, samples)
	outputFile.close()
	outputFile = open('mcmc/'+dir+'/sampler_chain.dat','w')
	numpy.save(outputFile, sampler.chain)
	outputFile.close()
	outputFile = open('mcmc/'+dir+'/sampler_lnprobability.dat','w')
	numpy.save(outputFile, sampler.lnprobability)
	outputFile.close()

	fp = open('mcmc/'+dir+'/chi.dat','w')
	probs=sampler.lnprobability.reshape([-1])
	for i,prob in enumerate(probs):
		print >>fp,i,-2.0*prob," ".join(str(x) for x in samples[i,:])
	fp.close()

	os.system('cp mcee.py mcmc/'+dir+'/')


def set_params(x,name):
	data="""output	0
source	1822
Lmin	3e31
Lscale	0.15
Tc	2e7
mass	1.2
radius	12.0
Bfield 	1e14
Edep	3.0
rhob	1e14
rhot	1e7
ytop	1e10
angle_mu	1
envelope	1
precalc	0
ngrid	50
SFgap	1
piecewise	0
timetorun	2000.0
Qimp	1
"""
	# Q,Lscale,Edep,Tc,M,R
	Lmin = PYL(x[3]*1e7, 1e14, x[4], x[5])
	
	params = {
		'Tc':x[3]*1e7,
		'Qimp':10.0**x[0],
		'Lscale':x[1],
		'Edep':x[2],
#		'Tt':x[2]*1e8,
#		'mdot':x[3],
		'mass':x[4],
		'radius':x[5],
		'Lmin':Lmin
#		'extra_Q':x[0],
#		'extra_y':x[1]
	}
	
	for key,value in params.items():
		data = re.sub("%s(\\t?\\w?).*\\n" % key,"%s\\t%g\\n" % (key,value),data)

	fout = open('/tmp/init.dat.'+name,'w')
	fout.write(data)
	fout.close()


def get_chisq(x):
	name = str(uuid.uuid4())
	set_params(x,name)
	# give crustcool a second parameter so that it looks in /tmp for the init.dat file
	data = subprocess.check_output( ["crustcool",name,'1'])
	chisq = float(re.search('chisq = ([-+]?[0-9]*\.?[0-9]+)',data).group(1))
#	print x[0],x[1],x[2],x[3],x[4],x[5],x[6],chisq
	os.system('rm /tmp/init.dat.'+name)
	return chisq


def lnprob(x):
	# minimum and maxiumum allowed values
	# (assume a flat prior within this range)

	# Q,Lscale,Edep,Tc,M,R
	xmin=numpy.array([-2.0,0.01,0.1,1.0,1.1,8.0])
	xmax=numpy.array([2.0,0.5,100.0,3.0,2.4,16.0])

#	xmin=numpy.array([1.0,-3.0,0.1,0.0,1.1,8.0])
#	xmax=numpy.array([100.0,3.0,100.0,3.0,2.5,16.0])
	#xmin=numpy.array([0.0,-3.0,0.0])
	#xmax=numpy.array([100.0,3.0,100.0])
	#xmin=numpy.array([0.0,-3.0,0.0,0.0])
	#xmax=numpy.array([100.0,3.0,100.0,1.0])
	if (len((x<xmin).nonzero()[0])>0 or len((x>xmax).nonzero()[0])>0):
		return -numpy.inf
	# causal limit  radius has to be R > 4.36km (M/Msun)
	if (x[5] < 4.36*x[4]):
		return -numpy.inf
	#print 'Trying ',x
	chisq=get_chisq(x)
	#print x[0],x[1],x[2],chisq
	return -chisq/2.0



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
