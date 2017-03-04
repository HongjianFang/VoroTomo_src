nthread = 20 
import numpy as np
import os
for ii in range(nthread):
	if not os.path.isdir(''.join(['fmm3drun',str(ii)])):
        	os.mkdir(''.join(['fmm3drun',str(ii)]))

with open('sources.in','r') as src:
	nsrc = int(src.readline())
	nsrcp = nsrc/nthread*np.ones((nthread,))
	nsrcp[-1] = nsrc-nsrc/nthread*(nthread-1)
	for ii in range(nthread):
	  with open(''.join(['fmm3drun',str(ii),'/sources.in']),'w') as srcp:
	   srcp.write(' '.join([str(int(nsrcp[ii])),'\n']))
	   for jj in range(int(nsrcp[ii])):
		line = src.readline()
		srcp.write(line)
		line = src.readline()
		srcp.write(line)
		npath = int(src.readline())
		srcp.write(' '.join([str(npath),'\n']))
		for kk in range(3*npath):
	 	 	line = src.readline()
			srcp.write(line)

aver = float(nsrc/nthread)
nrcp = np.zeros((nthread,),dtype=int)
with open('otimes.dat','r') as otime:
	nd = int(otime.readline())		  
	for ii in range(nd-int(aver)):
	  srcno = int(otime.readline().split()[1])
	  nrcp[int(np.ceil(srcno/aver))-1] = nrcp[int(np.ceil(srcno/aver))-1]+1
	  if nrcp[-1] == 1:
	    break 

nrcp[-1] = nd-np.sum(nrcp[:-1])	

with open('receivers.in','r') as rc:
	rc.readline()
	for ii in range(nthread):
	  with open(''.join(['fmm3drun',str(ii),'/receivers.in']),'w') as rcp:
	    rcp.write(' '.join([str(int(nrcp[ii])),'\n']))
	    for jj in range(int(nrcp[ii])):
	      line = rc.readline()
	      rcp.write(line)
	      line = rc.readline()
	      rcp.write(line)
	      line = int(rc.readline())
	      rcp.write(' '.join(['     '*2,str(line-int(np.sum(nsrcp[0:ii]))),'\n']))
	      line = rc.readline()
	      rcp.write(line)

#from shutil import copyfile
for ii in range(nthread):
	os.system(' '.join(['cp','frechgen.in','vgridsref.in','interfaces.in','mode_set.in','propgrid.in','vgrids.in','invert3d.in',''.join(['fmm3drun',str(ii)])]))
	os.system(' '.join(['cp',''.join(['fmm3drun',str(ii),'/sources.in']),''.join(['fmm3drun',str(ii),'/sourcesref.in'])]))

import threading
class mythread (threading.Thread):
	def __init__(self,threadID):
          threading.Thread.__init__(self)
	  self.threadID = threadID
	def run(self):
	  #print self.threadID
	  cmd1 = ''.join(['cd fmm3drun',str(self.threadID)])
	  cmd2 = 'frechgen'
	  cmd3 = ' '.join(['fm3d',''.join(['>fmm3dlog.out',str(self.threadID)])])
	  os.system(';'.join([cmd1,cmd2,cmd3]))
	  #print 'I am finished from', self.threadID

threads = []
for ii in range(nthread):
  t = mythread(ii)
  t.start()
  threads.append(t)

for t in threads:
  t.join()

print 'exiting all fmm3d'

with open('vgrids.in','r') as vgrids:
  vgrids.readline()
  nx,ny,nz =[int(ii) for ii in vgrids.readline().split()]
ngrid = nx*ny*nz*2
arrival = open('arrivals.dat','w')
nzero = 0
nrow = 0
from netCDF4 import Dataset
for ii in range(nthread):
	  with open(''.join(['fmm3drun',str(ii),'/arrivals.dat']),'r') as arr:
 	   for ls in arr:
	     line = ls.split()
	     didx = int(line[0])+int(np.sum(nrcp[:ii]))
	     sidx = int(line[1])+int(np.sum(nsrcp[:ii]))
	     line[0] = str(didx)
	     line[1] = str(sidx)
	     line.append('\n')
	     arrival.write('    '.join(line))
	  ds = Dataset(''.join(['fmm3drun',str(ii),'/frechet.nc']))
	  nrow = nrow + ds.dimensions.values()[0].size
	  nzero = nzero + ds.dimensions.values()[1].size
arrival.close()	      
ds.close()
nrows = np.zeros((nrow,),dtype=int)
nzeros = np.zeros((nzero,),dtype=float)
nzeros_id = np.zeros((nzero,),dtype=int)
ncount = 0
ncount2 = 0
for ii in range(nthread):
	  ds = Dataset(''.join(['fmm3drun',str(ii),'/frechet.nc']))
	  nrow1 = ds.dimensions.values()[0].size 
	  nrowseg = np.zeros((nrow1,),dtype=int)
	  nrowseg = ds.variables['Non_row'][:]
	  nzero1 =ds.dimensions.values()[1].size 
	  nzeroidseg = np.zeros((nzero1,),dtype=int)
	  nzeroidseg = ds.variables['Non_id'][:]
	  nzeroseg = np.zeros((nzero1,),dtype=float)
	  nzeroseg = ds.variables['Non_value'][:]
	  ncount3 = 0
	  for jj in range(nrow1):
		ncount = ncount+1
	        nrows[ncount-1] = nrowseg[jj] 
	        for kk in range(nrows[ncount-1]):
		   ncount2 = ncount2+1
		   ncount3 = ncount3+1
 		   nzeros[ncount2-1] = nzeroseg[ncount3-1]
		   tmp =  nzeroidseg[ncount3-1]
		   if tmp < ngrid:
		     nzeros_id[ncount2-1] = tmp
		   else:
		     nzeros_id[ncount2-1] = tmp+4*int(np.sum(nsrcp[:ii]))
	  ds.close()
#	  print 'finished',ii

dsnew = Dataset('frechet.nc','w',format='NETCDF3_CLASSIC')
nrowdim = dsnew.createDimension('nrow',nrow)
nzerodim = dsnew.createDimension('nzero',nzero)
ncnrow = dsnew.createVariable('Non_row',np.int32,('nrow',))
ncnzero = dsnew.createVariable('Non_value',np.float32,('nzero',))
ncnzero_id = dsnew.createVariable('Non_id',np.int32,('nzero',))
ncnrow[:] = nrows
ncnzero[:] = nzeros
ncnzero_id[:] = nzeros_id
dsnew.close()
