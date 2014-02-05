import numpy
import matplotlib
matplotlib.use('macosx')
import matplotlib.pyplot as plt
import triangle


fp = open('samples.dat','r')
samples = numpy.load(fp)
fp.close()


print samples.shape

samples = samples[:10000,:]

print samples.shape

fig = triangle.corner(samples,labels=[r"$T_{c,7}$", r"$Q_{imp}$", r"$T_{b,8}$",
						r"$\dot M$", r"$M (M_\odot)$", r"$R (km)$", r"$k_{n,crit}$"],
		quantiles=[0.16, 0.5, 0.84], plot_datapoints=True,bins=25,plot_ellipse=True)

fig.savefig('mcplot3.png')
