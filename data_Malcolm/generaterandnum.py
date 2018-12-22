import random
import numpy as np

with open('otimes.dat') as f:
	ndata = int(f.readline())
	eqidx = np.zeros((ndata,2))
	for ii in xrange(ndata):
		ln = f.readline().split()
		eqidx[ii,0] = int(ln[6])
		eqidx[ii,1] = int(ln[2])

nqall = int(np.max(eqidx[:,0])+1)
#uniqev,dataev = np.unique(eqdix[:,0],return_counts=True)
nchoose = 1000 
nnets = 100
choose = np.zeros(nchoose,dtype='int')
#dataidx = np.arange(ndata)+1
for ii in range(1,nnets+1):
  choose = random.sample(xrange(nqall),nchoose)
  #choose = np.sort(choose)
  choosenew = []
  eqidxnew = []
  phasenew = []
  for jj in range(nchoose):
    tmp = np.where(eqidx[:,0] == choose[jj])[0]
    choosenew = np.hstack([choosenew,tmp]) 
  finalidx = np.array(np.sort(choosenew),int)
  choosenew = finalidx+1
  eqidxnew = eqidx[finalidx,0]
  phasenew = eqidx[finalidx,1]

  choosenew = np.hstack([len(choosenew),choosenew])
  phasenew = np.hstack([len(choosenew),phasenew])
  tmp,eqidxnew_idx = np.unique(np.array(eqidxnew,dtype=int),return_inverse=True)
  eqidxnew = np.sort(np.unique(np.array(eqidxnew,dtype=int)))
  print ii,len(np.unique(eqidxnew)),len(np.unique(eqidxnew_idx))
  eqidxnew_idx = np.hstack([len(choosenew),eqidxnew_idx])
  choosenew = np.vstack([choosenew,eqidxnew_idx,phasenew])
  np.savetxt('tempdata/randomchoose'+str(ii)+'.dat',choosenew.T,fmt='%d')
  eqidxnew = np.hstack([len(eqidxnew),eqidxnew])
  np.savetxt('tempdata/updatesrc'+str(ii)+'.dat',eqidxnew,fmt='%d')
