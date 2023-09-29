import numpy
import matplotlib
# matplotlib.use('macosx')
import matplotlib.pyplot as plt
import corner
import sys
import numpy
import ns

def main():

	# get the directory name from the command line args
	dir = ''
	if len(sys.argv) > 1:
		dir += sys.argv[1]
		dir += '/'

	samples=read_samples(dir)
	if 1:
		plot_triangle(samples,dir)

	if 0:
		# 1D histograms of gravity
		samples=read_samples('1822_5/')
		#plot_hist1d(samples,'r','1822_run1',0)
		samples=read_samples('1822_6/')
		#plot_hist1d(samples,'b','1822_run2',0)
		plt.xlim([0.5,3.0])
		plt.ylim([0,0.2])
		samples=read_samples('1822_7/')
		#plot_hist1d(samples,'g','1822_run3',0)
		samples=read_samples('1822_8/')
		#plot_hist1d(samples,'r','1822_run4',0)
		samples=read_samples('1822_9/')
		#plot_hist1d(samples,'m','1822_run5',1)
		samples=read_samples('1822_10/')
		#plot_hist1d(samples,'c','1822_run6',1)
		samples=read_samples('1822_11/')
		plot_hist1d(samples,'k','1822_run7',1)

		samples=read_samples('1822_12/')
		plot_hist1d(samples,'r','1822_run8',1)

		samples=read_samples('1822_13/')
		plot_hist1d(samples,'b','1822_run9',1)

		samples=read_samples('1822_14/')
		plot_hist1d(samples,'g','1822_run10',1)

			
		plt.xlabel(r'$g_{14}$')
		plt.legend()
		plt.savefig('mcmc/'+dir+'g14_1822.pdf')

	if 0:
		# 1D histograms of gravity
		samples=read_samples('13/')
		plot_hist1d(samples,'r','T=54 eV')
		plt.xlim([0.5,3.0])
		plt.ylim([0,0.13])
		samples=read_samples('14/')
		plot_hist1d(samples,'b','T=63.1 eV')
		plt.xlabel(r'$g_{14}$')
		samples=read_samples('12/')
		plot_hist1d(samples,'k','current data')
		plt.legend()
		plt.savefig('g14.pdf')

	if 0:
		# 1D histograms of gravity
		samples=read_samples('13/')
		plot_hist1d(samples,'r','T=54 eV')
		plt.xlim([8.0,16.0])
		plt.ylim([0,0.07])
		samples=read_samples('14/')
		plot_hist1d(samples,'b','T=63.1 eV')
		plt.xlabel('R (km)')
		samples=read_samples('12/')
		plot_hist1d(samples,'k','current data')
		plt.legend(loc=2)
		plt.savefig('R.pdf')

	if 0:
		# 1D histograms of gravity
		samples=read_samples('13/')
		plot_hist1d(samples,'r','T=54 eV')
		plt.xlim([1.1,2.4])
		plt.ylim([0,0.14])
		samples=read_samples('14/')
		plot_hist1d(samples,'b','T=63.1 eV')
		plt.xlabel('R (km)')
		samples=read_samples('12/')
		plot_hist1d(samples,'k','current data')
		plt.legend()
		plt.savefig('mass.pdf')

	if 0:
		# 2D contours for Tc and gravity
		samples=read_samples('1822_5/')
		plot_hist(samples,'r')

		plt.savefig('g14_Q_1822.pdf')

	if 0:
		# 2D contours for Tc and gravity
		plt.subplot(3,3,1)
		
		samples=read_samples('1822_5/')
		plot_hist(samples,'r',0)
		
		plt.subplot(3,3,2)
		samples=read_samples('1822_6/')
		plot_hist(samples,'r',0)

		plt.subplot(3,3,3)
		samples=read_samples('1822_7/')
		plot_hist(samples,'r',0)		
		
		plt.subplot(3,3,4)
		samples=read_samples('1822_8/')
		plot_hist(samples,'r',0)		
		
		plt.subplot(3,3,5)
		samples=read_samples('1822_9/')
		plot_hist(samples,'r',1)		
		
		plt.subplot(3,3,6)
		samples=read_samples('1822_10/')
		plot_hist(samples,'r',1)		

		plt.subplot(3,3,7)
		samples=read_samples('1822_11/')
		plot_hist(samples,'r',1)		

		plt.subplot(3,3,8)
		samples=read_samples('1822_12/')
		plot_hist(samples,'r',1)		
		

		
		plt.savefig('MR_1822.pdf')




	if 0:
		# 2D contours for Tc and gravity
		plt.subplot(2,2,1)
		samples=read_samples('12/')
		plot_hist(samples,'r')

		plt.subplot(2,2,2)
		samples=read_samples('13/')
		plot_hist(samples,'r')

		plt.subplot(2,2,4)
		samples=read_samples('14/')
		plot_hist(samples,'r')

		plt.savefig('g14_Q_1822.pdf')


def read_samples(dir):
	
	# read in the samples
	fp = open('mcmc/'+dir+'samples.dat','rb')
	samples = numpy.load(fp)
	fp.close()
	# find out how many parameters
	nparams = samples.shape[1]
	# number of points
	n = samples.shape[0]
	# remove burn in points
	samples = samples[:int(1*n/10.0),:]
	
	return samples


def plot_hist1d(samples,col_string,label_string,flag):

	n = samples.shape[0]
	M = samples[:,4]
	if flag==0:
		R = samples[:,5]
		g = (6.67e-8*2e33/(1e14*1e10)) * M / R**2
		zz = 1.0/(1.0 - (2.0*1e14*g*1e5*R/9e20))**0.5
		g*=zz
	else:
		g=samples[:,5]
	#M = numpy.random.rand(n)*(2.5-1.1) + 1.1
	#R = numpy.random.rand(n)*(16.0-8.0) + 8.0
	Q = samples[:,1]
	mdot = samples[:,3]
	Tc = samples[:,0]

	if 0:
		ind = (mdot>0.05).nonzero()
		T_c = Tc[ind]
		g = g[ind]
		Q = Q[ind]
		R = R[ind]
		M = M[ind]
		mdot = mdot[ind]
		n=len(mdot)

#	n1, bins, patches = plt.hist(M, 50, weights=[1.0/n]*n,normed=False,alpha=1.0,color=col_string,histtype='step',label=label_string)
#	n1, bins, patches = plt.hist(R, 50, weights=[1.0/n]*n,normed=False,alpha=1.0,color=col_string,histtype='step',label=label_string)
	n1, bins, patches = plt.hist(g, 50, weights=[1.0/n]*n,normed=False,alpha=1.0,color=col_string,histtype='step',label=label_string)

def plot_hist(samples,col_string,flag):
	
	M = samples[:,4]
	if flag==0:
		R = samples[:,5]
	else:
		R=[]
		g=samples[:,5]
		for i,grav in enumerate(g):
			R.append(ns.R(M[i],grav))
		R=numpy.array(R)
	
	#Q = samples[:,1]
	#mdot = samples[:,3]
	#Tc = samples[:,0]

	Q = samples[:,0]
	Tc = samples[:,3]
	
	#n, bins, patches = plt.hist(g, 25, normed=1,alpha=0.1,color=col_string)
	
	if 0:
		ind = (mdot>0.05).nonzero()
		T_c = Tc[ind]
		g = g[ind]
		Q = Q[ind]
		mdot = mdot[ind]
	
	if 0:
		corner.hist2d(Tc,g,plot_datapoints=False,bins=25,extent=[(1,3),(0.5,3)])
		#triangle.hist2d(Tc,g,plot_datapoints=False,bins=25,extent=[(2,12),(0.5,3)])
		plt.xlabel(r'$T_c (10^7 K)$')
		plt.ylabel(r'$g_{14}$')

	if 0:
		corner.hist2d(g,Q,plot_datapoints=False,bins=25,extent=[(0.5,3),(-3,2)])
		plt.ylabel(r'$Q_{imp}$')
		plt.xlabel(r'$g_{14}$')

	if 1:
		corner.hist2d(R,M,plot_datapoints=False,bins=25,extent=[(8,16),(1.1,2.5)])
		plt.ylabel(r'$M$')
		plt.xlabel(r'$R$')

		Ranal = numpy.arange(100)*0.01*12.0 + 4.0
		Manal = []
		for R in Ranal:
			M=0.0
			while (ns.grav(M,R)<2.0):
				M+=0.01
			Manal.append(M)
		plt.plot(Ranal,Manal)
		Manal = []
		for R in Ranal:
			M=0.0
			while (ns.grav(M,R)<1.5):
				M+=0.01
			Manal.append(M)
		plt.plot(Ranal,Manal)
		Manal = []
		for R in Ranal:
			M=0.0
			while (ns.grav(M,R)<1.0):
				M+=0.01
			Manal.append(M)
		plt.plot(Ranal,Manal)
		
		plt.plot(Ranal,Ranal/4.36,'.')
		

	if 0:
		corner.hist2d(g,mdot,plot_datapoints=False,bins=25,extent=[(0.5,3),(0.0,0.5)])
		plt.ylabel(r'$Accretion rate$')
		plt.xlabel(r'$g_{14}$')

	


def plot_triangle(samples,dir):

	nparams = samples.shape[1]

	if nparams == 3:
		fig = corner.corner(samples,labels=[r"$Q_{imp}$", r"$L_{scale}$",
						r"$E_{dep}$"], 
			quantiles=[0.16, 0.5, 0.84], plot_datapoints=True,bins=50,plot_ellipse=False)

	#if nparams == 4:
	#	fig = corner.corner(samples,labels=[r"$T_{c,7}$", r"$Q_{imp}$", r"$T_{b,8}$",
	#					r"$\dot M$"], extents=[(4,11),(-3,2),(2,6),(0,0.5)],
	#		quantiles=[0.16, 0.5, 0.84], plot_datapoints=False,bins=50,plot_ellipse=False)

	#if nparams == 3:
	#	fig = corner.corner(samples,labels=[r"$T_{c,7}$", r"$Q_{imp}$", r"$T_{b,8}$"],
	#		extents=[(4,11),(-3,2),(2,6)],
	#				quantiles=[0.16, 0.5, 0.84], plot_datapoints=False,bins=50,plot_ellipse=False)

	if nparams == 6:
		fig = corner.corner(samples,labels=[r"$Q_{imp}$", r"$L_{scale}$",
						r"$E_{dep}$", r"$T_{c,7}$", r"$M (M_\odot)$", r"$g_{14}$"],
				quantiles=[0.16, 0.5, 0.84], plot_datapoints=True,bins=25,plot_ellipse=False,
				)

	if nparams == 8:
		fig = corner.corner(samples,labels=[r"$Q_{imp}$", r"$L_{scale}$",
					r"$T_{outer}$", r"$T_{c,7}$", r"$M (M_\odot)$", r"$g_{14}$",r"$\rho_{heat}$",r"$T_{inner}$"],
			quantiles=[0.16, 0.5, 0.84], plot_datapoints=True,bins=25,plot_ellipse=False)


				
#	if nparams == 6:
#		fig = corner.corner(samples,labels=[r"$Q_{imp}$", r"$L_{scale}$",
#					r"$E_{dep}$", r"$T_{c,7}$", r"$M (M_\odot)$", r"$R (km)$"],
#				quantiles=[0.16, 0.5, 0.84], plot_datapoints=False,bins=25,plot_ellipse=False,
#								)
		
				
#				if nparams == 6:
#					fig = corner.corner(samples,labels=[r"$T_{c,7}$", r"$Q_{imp}$", r"$T_{b,8}$",
#							r"$\dot M$", r"$M (M_\odot)$", r"$R (km)$"],extents=[(3,11),(-3,2),(2,6),(0,0.5),(1.2,2.4),(8,16)],
#							quantiles=[0.16, 0.5, 0.84], plot_datapoints=False,bins=50,plot_ellipse=False)

	fig.savefig('mcmc/'+dir+'mcplot.png')


if __name__ == '__main__':
    main()
