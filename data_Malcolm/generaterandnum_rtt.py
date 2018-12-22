import random
import numpy as np

data = np.loadtxt('differential_P.txt')
data[:,0:3] -= 1
nchoose = 10000
nset = 100
for ii in range(nset):
	print ii,'set'
	chooseid = random.sample(xrange(len(data)),nchoose)
	with open('tempdata/randomchoose_rtt'+str(ii+1)+'.dat','w') as f:
		f.write(str(nchoose)+'\n')
		for jj in range(nchoose):
#			f.write(' '.join(str(e) for e in data[chooseid[jj]])+'\n')
			f.write('%d %d %d %7.3f %7.3f\n' % tuple(data[chooseid[jj]]))
