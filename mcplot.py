import numpy
import matplotlib
matplotlib.use('macosx')
import matplotlib.pyplot as plt
import triangle
import sys

def main():

	# get the directory name from the command line args
	dir = ''
	if len(sys.argv) > 1:
		dir += sys.argv[1]
		dir += '/'

	samples=read_samples(dir)
	plot_triangle(samples,dir)

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
		plt.savefig('g14_new.pdf')

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

		plt.savefig('Tc_g14_new.pdf')


def read_samples(dir):
	
	# read in the samples
	fp = open('mcmc/'+dir+'samples.dat','r')
	samples = numpy.load(fp)
	fp.close()
	# find out how many parameters
	nparams = samples.shape[1]
	# number of points
	n = samples.shape[0]
	# remove burn in points
	samples = samples[1.0*n/4:,:]
	
	return samples


def plot_hist1d(samples,col_string,label_string):

	n = samples.shape[0]
	M = samples[:,4]
	R = samples[:,5]
	#M = numpy.random.rand(n)*(2.5-1.1) + 1.1
	#R = numpy.random.rand(n)*(16.0-8.0) + 8.0
	g = (6.67e-8*2e33/(1e14*1e10)) * M / R**2
	zz = 1.0/(1.0 - (2.0*1e14*g*1e5*R/9e20))**0.5
	g*=zz
	Q = samples[:,1]
	mdot = samples[:,3]
	Tc = samples[:,0]

	n1, bins, patches = plt.hist(g, 50, weights=[1.0/n]*n,normed=False,alpha=1.0,color=col_string,histtype='step',label=label_string)

def plot_hist(samples,col_string):
	
	M = samples[:,4]
	R = samples[:,5]
	g = (6.67e-8*2e33/(1e14*1e10)) * M / R**2
	zz = 1.0/(1.0 - (2.0*1e14*g*1e5*R/9e20))**0.5
	g*=zz
	Q = samples[:,1]
	mdot = samples[:,3]
	Tc = samples[:,0]
	
	#n, bins, patches = plt.hist(g, 25, normed=1,alpha=0.1,color=col_string)
	
	
	triangle.hist2d(Tc,g,plot_datapoints=False,bins=25,extent=[(2,12),(0.5,3)])
	plt.xlabel(r'$T_c (10^7 K)$')
	plt.ylabel(r'$g_{14}$')
	#triangle.hist2d(M/R,Q,plot_datapoints=False,bins=25,extent=[(0,0.33),(-3,2)])
	#triangle.hist2d(g,Q,plot_datapoints=False,bins=30,extent=[(0.5,2.5),(-3,2)],
	#			color=col_string,plot_contours=False,plot_ellipse=True)
	


def plot_triangle(samples,dir):

	nparams = samples.shape[1]

	if nparams == 4:
		fig = triangle.corner(samples,labels=[r"$T_{c,7}$", r"$Q_{imp}$", r"$T_{b,8}$",
						r"$\dot M$"], extents=[(4,11),(-3,2),(2,6),(0,0.5)],
			quantiles=[0.16, 0.5, 0.84], plot_datapoints=False,bins=50,plot_ellipse=False)

	if nparams == 3:
		fig = triangle.corner(samples,labels=[r"$T_{c,7}$", r"$Q_{imp}$", r"$T_{b,8}$"],
			extents=[(4,11),(-3,2),(2,6)],
					quantiles=[0.16, 0.5, 0.84], plot_datapoints=False,bins=50,plot_ellipse=False)

	if nparams == 6:
		fig = triangle.corner(samples,labels=[r"$T_{c,7}$", r"$Q_{imp}$", r"$T_{b,8}$",
				r"$\dot M$", r"$M (M_\odot)$", r"$R (km)$"],extents=[(3,11),(-3,2),(2,6),(0,0.5),(1.2,2.4),(8,16)],
				quantiles=[0.16, 0.5, 0.84], plot_datapoints=False,bins=50,plot_ellipse=False)

	fig.savefig('mcmc/'+dir+'mcplot.png')


if __name__ == '__main__':
    main()
