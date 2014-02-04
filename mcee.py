import os
import random
import numpy
import matplotlib
matplotlib.use('macosx')
import matplotlib.pyplot as plt
import emcee
import triangle
import uuid
import time


def set_params(x,name):
	fin = open('init/init.dat.mcmc','r')
	fout = open('init/init.dat.'+name,'w')
	data_in = fin.read()
	
	params = {
		'Tc':x[0]*1e7,
		'Qimp':10.0**x[1],
		'Tt':x[2]*1e8,
		'mdot':x[3],
		'mass':x[4],
		'radius':x[5],
		'gpe':0,
		'toutburst':12.5
	}
	
	for line in data_in.split('\n'):
		flag = 0
		for key,value in params.items():
			if line.find(key) != -1:
				fout.write("%s\t%g\n" % (key,value))
				flag = 1
		if flag==0:
			fout.write(line+'\n')
	fin.close()
	fout.close()


def get_chisq(x):
	name = str(uuid.uuid4())
	set_params(x,name)
	os.system('crustcool '+name+' >tmp'+name)
	f=open('tmp'+name,'r')
	data = f.read()
	for line in data.split('\n'):
		if line.find('chisq') != -1 and line.find('chisq_nu') == -1:
			chisq = float(line.split('=')[1])
	f.close()
#	print x[0],x[1],x[2],chisq
#	print x[0],x[1],x[2],x[3],x[4],x[5],chisq
	os.system('rm tmp'+name)
	os.system('rm init/init.dat.'+name)	
	return chisq


def lnprob(x):
	if (x[4] < 1.1 or x[4]>2.5 or x[5]<8.0 or x[5]>16.0 or x[3]<0.0 or x[3]>1.0):
		return -numpy.inf
	if (x[1]<-3.0 or x[1]>2.0 or x[2]<0.0 or x[0]<0.0):
		return -numpy.inf
	return -get_chisq(x)/2.0



nwalkers=200
ndim=6
# parameters are   x = [Tc7, Qimp, Tb8, mdot, M, R]
p0 = emcee.utils.sample_ball([7.0,0.0,5.0,0.1,1.62,11.2],[0.1,0.1,0.1,0.01,0.1,0.5],nwalkers)
#p0 = emcee.utils.sample_ball([3.5,0.0,4.2,0.1,1.62,11.2],[0.1,0.1,0.1,0.01,0.1,0.5],nwalkers)
#p0 = emcee.utils.sample_ball([3.5,0.0,4.2],[0.1,0.1,0.1],nwalkers)

sampler=emcee.EnsembleSampler(nwalkers,ndim,lnprob,threads=15)

start_time = time.time()
sampler.run_mcmc(p0, 120)
print 'time to run = ',time.time() - start_time,'seconds'

# throw away the first nburn steps as burn-in
nburn = 20
samples=sampler.chain[:,nburn:,:].reshape((-1,ndim))
fig = triangle.corner(samples,labels=[r"$T_{c,7}$", r"$Q_{imp}$", r"$T_{b,8}$",
						r"$\dot M$", r"$M (M_\odot)$", r"$R (km)$"],
		quantiles=[0.16, 0.5, 0.84])
#fig = triangle.corner(samples,labels=[r"$T_{c,7}$", r"$Q_{imp}$", r"$T_{b,8}$"],
#		 quantiles=[0.16, 0.5, 0.84])
fig.savefig("triangle.pdf")

print("Mean acceptance fraction: {0:.3f}"
                .format(numpy.mean(sampler.acceptance_fraction)))

outputFile = open('samples.dat','w')
numpy.save(outputFile, samples)
outputFile.close()
