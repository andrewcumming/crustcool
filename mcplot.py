import numpy
import matplotlib
matplotlib.use('macosx')
import matplotlib.pyplot as plt
import triangle
import sys

dir = ''
if len(sys.argv) > 1:
	dir += sys.argv[1]
	dir += '/'

fp = open('mcmc/'+dir+'samples.dat','r')
samples = numpy.load(fp)
fp.close()

nparams = samples.shape[1]
n = samples.shape[0]
samples = samples[1.0*n/4:,:]

#fig = triangle.corner(samples,labels=[r"$T_{c,7}$", r"$Q_{imp}$", r"$T_{b,8}$",
#						r"$\dot M$", r"$M (M_\odot)$", r"$R (km)$", r"$k_{n,crit}$"],
#		quantiles=[0.16, 0.5, 0.84], plot_datapoints=True,bins=24,plot_ellipse=False)
#fig = triangle.corner(samples,labels=[r"$Q\ (MeV)$", r"$y\ (g\ cm^{-2})$"],
#		quantiles=[0.16, 0.5, 0.84], plot_datapoints=True,bins=24,plot_ellipse=False)

if nparams == 4:
	fig = triangle.corner(samples,labels=[r"$T_{c,7}$", r"$Q_{imp}$", r"$T_{b,8}$",
					r"$\dot M$"], extents=[(4,11),(-3,2),(2,6),(0,0.5)],
		quantiles=[0.16, 0.5, 0.84], plot_datapoints=False,bins=50,plot_ellipse=False)

if nparams == 3:
	fig = triangle.corner(samples,labels=[r"$T_{c,7}$", r"$Q_{imp}$", r"$T_{b,8}$"],
		extents=[(4,11),(-3,2),(2,6)],
				quantiles=[0.16, 0.5, 0.84], plot_datapoints=False,bins=50,plot_ellipse=False)

fig.savefig('mcmc/'+dir+'mcplot.png')
