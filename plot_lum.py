import sys
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
sns.set_style("ticks",{'xtick.direction': 'in','ytick.direction': 'in'})
sns.set_context("paper",font_scale=1.6)
fig = plt.figure(figsize=(8,6))
ax = fig.add_subplot(1,1,1)

def read_lightcurve(name):
	f = open("out/prof"+name)
	t=np.array([])
	L=np.array([])
	for line in f:
		# convert to time in days
		t=np.append(t,float(line.split()[0])/(24.0*3600.0))
		# convert temperature to eV
		L=np.append(L,float(line.split()[1]))

		# Note that the output from the code is already redshifted,
		# we don't need to do it here
	f.close()
	return t,L

# The source name can be given on the command line
if len(sys.argv)>1:
	source_name = str(sys.argv[1])
else:
	source_name = '1659'

if 1:
	area = 4.0*np.pi*1.2e6**2

	#t,L = read_lightcurve("_1731_He")
	#plt.plot(t,L*area,'k:',label=r'$\mathregular{He\ envelope}$')
	
	#t,L = read_lightcurve("_1731_Fe")
	#plt.plot(t,L*area,'k--',label=r'$\mathregular{Fe\ envelope}$')

	t,L = read_lightcurve("")
	plt.plot(t,L*area,'r')

	plt.legend()

# Axis labels etc. and output
plt.ylabel('$L \mathrm{(erg\ s^{-1})}$')
plt.xlabel('Time after outburst (days)')
plt.xlim((1.0,1e4))
plt.ylim((1e31,1e35))
#plt.ylim((58,110))
ax.set_xscale('log')
ax.set_yscale('log')
#plt.savefig("lum.pdf")
plt.show()
