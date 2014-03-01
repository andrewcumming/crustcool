import random
import numpy
import matplotlib
matplotlib.use('macosx')
import matplotlib.pyplot as plt
import emcee
import triangle
import uuid
import time
import subprocess
import re
import os
import datetime

def set_params(x,name):
	data="""output	0
mass	1.2
radius	13.0
Bfield 	0
mdot	1.0
precalc	1
ngrid	50
SFgap	1
kncrit	0.0
sph	0
piecewise	0
timetorun	1000.0
neutrinos	1
instant	0
toutburst	1.5
accreted	1
gpe	0
ytop	1e10
Tt	4.2e8
Qimp	1
Tc	1e8
cooling_bc	1
extra_heating	1
extra_Q	3.0
extra_y	13.0
"""
	
	params = {
#		'Tc':x[0]*1e7,
#		'Qimp':10.0**x[1],
#		'Tt':x[2]*1e8,
#		'mdot':x[3],
#		'mass':x[4],
#		'radius':x[5]
		'extra_Q':x[0],
		'extra_y':x[1]
	}
	
	for key,value in params.items():
		data = re.sub("%s(\\t?\\w?).*\\n" % key,"%s\\t%f\\n" % (key,value),data)

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
	#xmin=numpy.array([0.0,-3.0,0.0,0.0,1.1,8.0])
	#xmax=numpy.array([100.0,3.0,100.0,3.0,2.5,16.0])
	xmin=numpy.array([0.0,11.0])
	xmax=numpy.array([10.0,15.0])
	if (len((x<xmin).nonzero()[0])>0 or len((x>xmax).nonzero()[0])>0):
		return -numpy.inf
	chisq=get_chisq(x)
	#print x[0],x[1],chisq
	return -chisq/2.0


nwalkers, ndim = 100, 2

# parameters are   x = [Tc7, Qimp, Tb8, mdot, M, R, kncrit]
p0 = emcee.utils.sample_ball([3.0,13.0],[1.0,1.0],nwalkers)
#p0 = emcee.utils.sample_ball([20.0,1.0,1.5,2.0,1.2,13.0],[0.3,0.1,0.3,0.2,0.01,0.5],nwalkers)

sampler=emcee.EnsembleSampler(nwalkers,ndim,lnprob,threads=15)

print 'Starting run at ', str(datetime.datetime.now())
start_time = time.time()
sampler.run_mcmc(p0, 100)
print 'time to run = ',time.time() - start_time,'seconds'
print("Mean acceptance fraction: {0:.3f}"
                .format(numpy.mean(sampler.acceptance_fraction)))

# throw away the first nburn steps as burn-in
nburn = 20
samples=sampler.chain[:,nburn:,:].reshape((-1,ndim))
fig = triangle.corner(samples)
#		fig = triangle.corner(samples,labels=[r"$T_{c,7}$", r"$Q_{imp}$", r"$T_{b,8}$",
#								r"$\dot M$", r"$M (M_\odot)$", r"$R (km)$", r"k_{n,crit}"],
#				quantiles=[0.16, 0.5, 0.84])
fig.savefig("mcmc/triangle.pdf")

outputFile = open('mcmc/samples.dat','w')
numpy.save(outputFile, samples)
outputFile.close()
outputFile = open('mcmc/sampler_chain.dat','w')
numpy.save(outputFile, sampler.chain)
outputFile.close()
outputFile = open('mcmc/sampler_lnprobability.dat','w')
numpy.save(outputFile, sampler.lnprobability)
outputFile.close()
