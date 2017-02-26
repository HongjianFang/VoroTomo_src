nthread = 16 
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
nrcp = np.zeros(nthread,)
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
frech = open('frechet.dat','w')
arrival = open('arrivals.dat','w')
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
	  with open(''.join(['fmm3drun',str(ii),'/frechet.dat']),'r') as fre:
	   for kk in range(int(nrcp[ii])):
	     line = fre.readline().split()
	     didx = int(line[0])+int(np.sum(nrcp[:ii]))
	     sidx = int(line[1])+int(np.sum(nsrcp[:ii]))
	     line[0] = str(didx)
	     line[1] = str(sidx)
	     line.append('\n')
	     frech.write('    '.join(line))
	     for jj in range(int(line[4])):
		line = fre.readline().split()
		if int(line[0])<=ngrid:
	          line.append('\n')
		  frech.write('  '.join(line))
		else:
		  spidx = int(line[0])+4*int(np.sum(nsrcp[:ii]))
		  line[0] = str(spidx)
	          line.append('\n')
		  frech.write('  '.join(line))
arrival.close()	      
frech.close()
