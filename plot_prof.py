import sys
import numpy
import matplotlib.pyplot as plt
import seaborn as sns
sns.set_style("ticks",{'xtick.direction': u'in','ytick.direction': u'in'})
sns.set_context("paper",font_scale=1.6)
fig = plt.figure(figsize=(8,6))
ax = fig.add_subplot(1,1,1)

def read_lightcurve(name):
	f = open("gon_out/out"+name)
	ngrid = int(f.readline().split()[0])
	f.readline()
	rho=numpy.array([])
	T=numpy.array([])
	for i in range(ngrid):
		line = f.readline()
		rho=numpy.append(rho,float(line.split()[5]))
		T=numpy.append(T,float(line.split()[1]))
	f.close()
	return rho,T

rho,T = read_lightcurve("")
plt.plot(rho,T,'r')


# Axis labels etc. and output
plt.xlabel(r'$\rho\ (\mathrm{g\ cm^{-3}})$')
plt.ylabel(r'$T\ (\mathrm{K})$')
#plt.xlim((30,1e4))
#plt.ylim((40,125))
#plt.ylim((58,110))
ax.set_xscale('log')
ax.set_yscale('log')
#plt.savefig("prof.pdf")
plt.show()
