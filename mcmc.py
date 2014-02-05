import os
import random
import numpy
import matplotlib
matplotlib.use('macosx')
import matplotlib.pyplot as plt
import pylab

pylab.ion()



def set_params(x,name):
	data="""output	0
mass	1.62
radius	11.2
Bfield 	0
mdot	0.1
precalc	1
ngrid	50
SFgap	5
kncrit	0.0
sph	0
piecewise	0
timetorun	4000.0
neutrinos	1
instant	0
toutburst	1.6
accreted	1
gpe	0
ytop	1e12
Tt	4.2e8
Qimp	3.5
Tc	3.0e7
"""
	
	params = {
		'Tc':x[0]*1e7,
		'Qimp':10.0**x[1],
		'Tt':x[2]*1e8,
	#	'mdot':x[3],
	#	'mass':x[4],
	#	'radius':x[5]
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

	chisq_new = get_chisq(anew)

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
		chisq = chisq_new
		
	Tc.append(a[0])
	Qimp.append(a[1])
	Tb.append(a[2])
	plt.plot(Tc[-1],Qimp[-1],'r+')
	plt.draw()
	fmc.write("%g %g %g %g\n" % (a[0],a[1],a[2],chisq))
	if count % 10 == 0:
		fmc.flush()
		
	print a[0],a[1], 10.0**a[1],a[2], chisq, accept_count, count, 1.0*accept_count/count

fmc.close()