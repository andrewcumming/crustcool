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
	fin = open('init/init.dat.mcmc','r')
	data = fin.read()
	fin.close()
	
	params = {
		'Tc':x[0]*1e7,
		'Qimp':10.0**x[1],
		'Tt':x[2]*1e8,
		'mdot':x[3],
		'mass':x[4],
		'radius':x[5],
		'gpe':0,
		'toutburst':1.6,
		'SFgap':5,
		'kncrit':0.0,
		'precalc':1
	}
	
	for key,value in params.items():
		data = re.sub("%s(\\t?\\w?).*\\n" % key,"%s\\t%f\\n" % (key,value),data)

	fout = open('init/init.dat.'+name,'w')
	fout.write(data)
	fout.close()


def get_chisq(x):
	name = str(uuid.uuid4())
	set_params(x,name)
	data = subprocess.check_output( ["crustcool",name])
	chisq = float(re.search('chisq = ([-+]?[0-9]*\.?[0-9]+)',data).group(1))
#	print x[0],x[1],x[2],chisq
#	print x[0],x[1],x[2],x[3],x[4],x[5],x[6],chisq
	os.system('rm init/init.dat.'+name)	
	return chisq


def lnprob(x):
	if (x[4] < 1.1 or x[4]>2.5 or x[5]<8.0 or x[5]>16.0 or x[3]<0.0 or x[3]>3.0):
		return -numpy.inf
	#if (x[6] < 0.0 or x[6] > 1.0):
	#	return -numpy.inf
	if (x[1]<-3.0 or x[1]>2.0 or x[2]<0.0 or x[0]<0.0):
		return -numpy.inf
	return -get_chisq(x)/2.0



nwalkers=200
ndim=6
# parameters are   x = [Tc7, Qimp, Tb8, mdot, M, R, kncrit]
#p0 = emcee.utils.sample_ball([3.5,0.0,4.2,0.1,1.62,11.2,0.5],[0.1,0.1,0.1,0.01,0.1,0.5,0.05],nwalkers)
p0 = emcee.utils.sample_ball([3.5,0.0,4.2,1.0,1.62,11.2],[0.3,0.1,0.3,0.2,0.1,0.5],nwalkers)
#p0 = emcee.utils.sample_ball([3.5,0.0,4.2,0.1,1.62,11.2],[0.1,0.1,0.1,0.01,0.1,0.5],nwalkers)
#p0 = emcee.utils.sample_ball([3.5,0.0,4.2],[0.1,0.1,0.1],nwalkers)

sampler=emcee.EnsembleSampler(nwalkers,ndim,lnprob,threads=15)

print 'Starting run at ', str(datetime.datetime.now())

start_time = time.time()
sampler.run_mcmc(p0, 220)
print 'time to run = ',time.time() - start_time,'seconds'

# throw away the first nburn steps as burn-in
nburn = 20
samples=sampler.chain[:,nburn:,:].reshape((-1,ndim))
fig = triangle.corner(samples,labels=[r"$T_{c,7}$", r"$Q_{imp}$", r"$T_{b,8}$",
						r"$\dot M$", r"$M (M_\odot)$", r"$R (km)$", r"k_{n,crit}"],
		quantiles=[0.16, 0.5, 0.84])
#fig = triangle.corner(samples,labels=[r"$T_{c,7}$", r"$Q_{imp}$", r"$T_{b,8}$"],
#		 quantiles=[0.16, 0.5, 0.84])
fig.savefig("triangle.pdf")

print("Mean acceptance fraction: {0:.3f}"
                .format(numpy.mean(sampler.acceptance_fraction)))

outputFile = open('samples.dat','w')
numpy.save(outputFile, samples)
outputFile.close()
