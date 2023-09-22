import numpy
import matplotlib.pyplot as plt
fig = plt.figure()
ax = fig.add_subplot(1,1,1)

def read_profile():
	timestamp = f.readline()
	if not timestamp:
		return timestamp,-1,-1
	rho=numpy.array([])
	T=numpy.array([])
	for i in range(ngrid):
		line = f.readline()
		rho=numpy.append(rho,float(line.split()[5]))
		T=numpy.append(T,float(line.split()[1]))
	timestamp = "%.2g" % (float(timestamp)/(24.0*3600.0))
	timestamp = 't = '+timestamp+' d'
	return timestamp,rho,T

output_png = False

f = open("out/out")
ngrid = int(f.readline().split()[0])

plt.ion()
# Axis labels etc. and output
plt.xlabel(r'$\rho\ (\mathrm{g\ cm^{-3}})$')
plt.ylabel(r'$T\ (\mathrm{K})$')
plt.xlim((3e8,2e14))
plt.ylim((1e7,5e8))
ax.set_xscale('log')
ax.set_yscale('log')

timestamp,rho,T = read_profile()
x1, = plt.plot(rho,T,'k')
ann = plt.annotate(timestamp, xy=(1e12,3e8))
fig.canvas.draw()

if output_png:
	i=1
	plt.savefig("png/%03d.png" % i)
while timestamp:
	timestamp,rho,T = read_profile()
	if timestamp:
		x1.set_ydata(T)
		ann.set_text(timestamp)
		fig.canvas.draw()
		if output_png:
			i = i + 1
			plt.savefig("png/%03d.png" % i)
		plt.pause(0.001)

eval(input('Press any key to end'))
f.close()
