import os
import random
import numpy
import matplotlib
matplotlib.use('macosx')
import matplotlib.pyplot as plt
import pylab

pylab.ion()

def set_params(a):
	fin = open('init/init.dat.mcmc','r')
	fout = open('init.dat','w')
	data_in = fin.read()
	
	Tc_to_set = a[0]*1e7
	Qimp_to_set = 10.0**a[1]
	Tt_to_set = a[2]*1e8
	
	for line in data_in.split('\n'):
		flag = 0
		if line.find('Tc') != -1:
			fout.write("Tc\t%g\n" % Tc_to_set)
			flag = 1
		if line.find('Qimp') != -1:
			fout.write("Qimp\t%g\n" % Qimp_to_set)
			flag = 1
		if line.find('Tt') != -1:
			fout.write("Tt\t%g\n" % Tt_to_set)
			flag = 1
		if flag==0:
			fout.write(line+'\n')
	fin.close()
	fout.close()


def get_chisq():
	os.system('crustcool >tmp')
	f=open('tmp','r')
	data = f.read()
	for line in data.split('\n'):
		if line.find('chisq') != -1 and line.find('chisq_nu') == -1:
			chisq = float(line.split('=')[1])
	f.close()
	return chisq


# a = [Tc7, Qimp, Tb8]
a=numpy.array([4.0,0.0,4.2])
set_params(a)
chisq = get_chisq()

count = 0
accept_count = 0

Tc=[]
Qimp=[]
Tb=[]

fmc=open('out/mcmc.dat','w')

fac = 0.3

while True:
	# make a jump
	anew = a.copy()
	# Tc  don't allow to go <1.0 (i.e. Tc<1e7K)
	anew[0]=a[0]+fac*0.2*random.gauss(0.0,1.0)
	while anew[0]<1.0:
		anew[0]=a[0]+fac*0.2*random.gauss(0.0,1.0)
	# Qimp:  only allow a range 1e-3 to 1000.0
	anew[1]=a[1]+fac*0.2*random.gauss(0.0,1.0)
	while anew[1]<-3.0 or anew[1]>3.0:
		anew[1]=a[1]+fac*0.2*random.gauss(0.0,1.0)
	# Ttop
	anew[2]=a[2]+fac*0.2*random.gauss(0.0,1.0)
	set_params(anew)
	chisq_new = get_chisq()

	accept_flag=0
	if chisq_new<chisq:
		accept_flag=1
	else:
		if random.random()<numpy.exp(-0.5*(chisq_new-chisq)):
			accept_flag=1

	count+=1
	if accept_flag:
		accept_count+=1
		a = anew
		Tc.append(a[0])
		Qimp.append(a[1])
		Tb.append(a[2])
		chisq = chisq_new
		plt.plot(Tc[-1],Qimp[-1],'r+')
		plt.draw()
		fmc.write("%g %g %g %g\n" % (a[0],a[1],a[2],chisq))
		if accept_count % 10 == 0:
			fmc.flush()
		
	print a[0],a[1], 10.0**a[1],a[2], chisq, accept_count, count, 1.0*accept_count/count

fmc.close()