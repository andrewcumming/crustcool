import sys
import numpy
import matplotlib.pyplot as plt
import seaborn as sns
sns.set_style("ticks",{'xtick.direction': 'in','ytick.direction': 'in'})
sns.set_context("paper",font_scale=1.6)
fig = plt.figure(figsize=(8,6))
ax = fig.add_subplot(1,1,1)

def read_lightcurve(name):
	f = open("out/prof"+name)
	t=numpy.array([])
	Teff=numpy.array([])
	for line in f:
		# convert to time in days
		t=numpy.append(t,float(line.split()[0])/(24.0*3600.0))
		# convert temperature to eV
		Teff=numpy.append(Teff,float(line.split()[4])*1.38e-16/1.6e-12)
		# Note that the output from the code is already redshifted,
		# we don't need to do it here
	f.close()
	return t,Teff

# The source name can be given on the command line
if len(sys.argv)>1:
	source_name = str(sys.argv[1])
else:
	source_name = '1659'

t,Teff = read_lightcurve("")
plt.plot(t,Teff,'r')

# Read in the data and plot it
f = open("data/"+source_name)
t0, n = f.readline().split()
t0 = float(t0)
t = numpy.array([])
Teff = numpy.array([])
Terr = numpy.array([])
for line in f:
	t=numpy.append(t,float(line.split()[0])-t0)
	Teff=numpy.append(Teff,float(line.split()[3]))
	Terr=numpy.append(Terr,float(line.split()[4]))
f.close()

plt.errorbar(t,Teff,fmt='ko',yerr=Terr,xerr=0.0,ecolor='k')

# Axis labels etc. and output
plt.ylabel('$\mathregular{T_{eff}^\infty\ (eV)}$')
plt.xlabel('Time after outburst (days)')
plt.xlim((30,1e4))
plt.ylim((40,125))
#plt.ylim((58,110))
ax.set_xscale('log')
plt.savefig("tc.pdf")
#plt.show()
