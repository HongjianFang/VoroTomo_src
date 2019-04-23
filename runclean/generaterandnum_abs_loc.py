def generate_abs_loc():
  import random
  import numpy as np
  
  with open('abs.dat') as f:
  	ndata = int(f.readline())
  	eqidxsrc = np.zeros(ndata,dtype=int)
  	eqidxrc = np.zeros(ndata,dtype=int)
  	eqidxabst = np.zeros(ndata,)
  	for ii in xrange(ndata):
  		ln = f.readline().split()
  		eqidxsrc[ii] = int(ln[0])-1
  		eqidxrc[ii] = int(ln[1])-1
  		eqidxabst[ii] = float(ln[2])
  		#eqidx[ii,2] = int(ln[2])
  
  with open('otimes.dat') as f:
	    ndataall = int(f.readline())
	    srcallidx = np.zeros(ndataall,dtype=int)
	    rcallidx = np.zeros(ndataall,dtype=int)
	    for ii in range(ndataall):
		    ln = f.readline().split()
		    srcallidx[ii] = int(ln[6])
		    rcallidx[ii] = int(ln[7])
  
  nchoose = int(ndata*0.2)
  choosenew = np.sort(random.sample(xrange(ndata),nchoose))
  with open('tempdata/randomchoosereloc_abs.dat','w') as f:
	    f.write(str(nchoose)+'\n')
   	    for jj in range(nchoose):
		ii = choosenew[jj]
	        idxsrc = np.where(rcallidx==eqidxrc[ii])[0]
	        idxrc = np.where(srcallidx[idxsrc]==eqidxsrc[ii])[0]
	        f.write('%d %d %7.4f \n' % (idxsrc[idxrc]+1,1,eqidxabst[ii]))
  print 'finishing random picking events'
