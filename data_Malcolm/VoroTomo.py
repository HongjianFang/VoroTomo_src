#!/usr/bin/env python
########################################################################################
# Script to run joint inversion with body wave and surface wave data
# Hongjian Fang @ USTC 2016.11.11
########################################################################################
#
import os
from shutil import copyfile,move
from netCDF4 import Dataset
from datetime import datetime
#import sys
import numpy as np
#

# be careful about the different paremeter to invert. the index might be wrong. Now this set is suit for Vp Vs and hypocent at the same time
# forward
nthread = 8 
fmm='../../bin/fm3d'
frech='frechgen'
ttim='arrivals.dat'
cvg='vgrids.in'
cig='interfaces.in'
csl='sources.in'
crc='receivers.in'
diagnos='fm3dlog.out'
# surface wave
#########################################################
surfdirect='surfdirect'
mtravsurf='mtimessurf.dat'
ttimessurf='ttimessurf.dat'
diagnossurf='surfdirect.out'
#
# inversion 
inv='../randmesh_src/relocation'
invinput = 'invert3d.in'
itn='inviter.in'
mtrav='mtimes.dat'
rtrav='rtimes.dat'
ivg='vgridsref.in'
irc='receiversori.in'
iig='interfacesref.in'
isl='sourcesref.in'
#
# residuals
resid='residuals'
resout='residuals.dat'
#
generatemesh='../randmesh_src/generatemesh'
subspaceproj='../randmesh_src/subspaceproj'
backproj_ind='../randmesh_src/backproj_ind'
combine_sub='../randmesh_src/combine_sub'
# outfiles
gmeshout='gmeshout.out'
backprojout='backproj.out'
combineout='combine.out'
subpout='subspacep.out'
ifile='tomo3d_rm.in'
#
if os.path.isfile(gmeshout):
	os.remove(gmeshout)
if os.path.isfile(backprojout):
	os.remove(backprojout)
if os.path.isfile(combineout):
	os.remove(combineout)
if os.path.isfile(subpout):
	os.remove(subpout)

import threading
class mythread (threading.Thread):
	def __init__(self,threadID):
          threading.Thread.__init__(self)
	  self.threadID = threadID
	def run(self):
	  #print self.threadID
	  cmd1 = ''.join(['cd fmm3drun',str(self.threadID)])
	  cmd2 = 'frechgen'
	  cmd3 = ' '.join([fmm,''.join(['>fmm3dlog.out',str(self.threadID)])])
	  os.system(';'.join([cmd1,cmd2,cmd3]))
	  #print 'I am finished from', self.threadID

class parallelrunsub (threading.Thread):
	def __init__(self,threadID):
          threading.Thread.__init__(self)
	  self.threadID = threadID
	def run(self):
	  inet = self.threadID
          cmd1 = subspaceproj+ ' '+str(inet)+' >> '+subpout
          cmd2 = backproj_ind+' '+str(inet)+' >> '+ backprojout
	  os.system(';'.join([cmd1,cmd2]))
	
def parallel_fm3d(invtype):
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
	srcidx = np.zeros(nd,)
  	for ii in range(nd):
 	  lns = otime.readline().split()
	  srcidx[ii] = int(lns[6])
  with open('otimes.dat','r') as otime:
  	nd = int(otime.readline())		  
  	for ii in range(nd-int(aver)):
 	  lns = otime.readline().split()
	  srcno = int(lns[1])
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
# note bug here, only for Vp and Vs together
  #ngrid = nx*ny*nz*2
# For Vp only
  ngrid = nx*ny*nz
  arrival = open('arrivals.dat','w')
  nzero = 0
  nrow = 0
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
  	  ds.close()
  arrival.close()	      
  
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
	  if invtype:
  	    for jj in range(nrow1):
  		ncount = ncount+1
  	        nrows[ncount-1] = nrowseg[jj] 
  	        for kk in range(nrows[ncount-1]):
  		   ncount2 = ncount2+1
  		   ncount3 = ncount3+1
   		   nzeros[ncount2-1] = nzeroseg[ncount3-1]
  		   tmp =  nzeroidseg[ncount3-1]
		   #if invtype:#nrows[ncount-1]>4 and tmp <= ngrid:
		   if tmp <= ngrid:
  		     nzeros_id[ncount2-1] = tmp
  		   else:
		     nzeros_id[ncount2-1] = ngrid+np.mod((tmp-1-ngrid),4)+1+4*srcidx[ncount-1]
	  else:
  	    for jj in range(nrow1):
  		ncount = ncount+1
  	        nrows[ncount-1] = nrowseg[jj] 
  	        for kk in range(nrows[ncount-1]):
  		   ncount2 = ncount2+1
  		   ncount3 = ncount3+1
   		   nzeros[ncount2-1] = nzeroseg[ncount3-1]
  		   tmp =  nzeroidseg[ncount3-1]
		   #if invtype:#nrows[ncount-1]>4 and tmp <= ngrid:
		   nzeros_id[ncount2-1] = np.mod((tmp-1),4)+1+4*srcidx[ncount-1]
  	  ds.close()
  if (ncount!=nd):
	print 'warning: something is wrong'
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
  return

##############################
def run_randmesh(runid,nomesh):
  with open(ifile, 'r') as fp:
      lines = fp.readlines()
  NI = int(lines[0].strip())
  BGITER = int(lines[1].strip())
  BINV = int(lines[2].strip())
  JOINT = int(lines[3].strip())
  numfiles = int(lines[4].strip())
  velinv = 1
  ###############################
  ###############################
  #
  # If necessary, copy the initial velocity,
  # interface and source files to the current files
  # read by fm3d. Generate the frechet.in file and
  # set the iteration number to 1.
  #
  start_time = datetime.now()
  if BGITER == 0:
    copyfile(ivg,cvg)
    copyfile(irc,crc)
    copyfile(iig,cig)
    copyfile(isl,csl)
    ITER=1
    with open(itn, 'w') as fp:
        fp.write(str(ITER))
    with open(invinput,'r') as fp:
        data = fp.readlines()
    data[22] = lines[ITER+4].split('|')[0].strip()+'\n'
    data[23] = lines[ITER+4].split('|')[1].strip()+'\n'
    data[24] = lines[ITER+4].split('|')[2].strip()+'\n'
    data[30] = lines[3]
    velinv = int(data[22].split()[0])
    with open(invinput,'w') as fp:
        fp.writelines(data)

    os.system(frech)

  
  #
  # Run FMM
  #
# uncomment this
  if BGITER==0 or (BGITER==1 and BINV==1 ):
    if(JOINT==0 or JOINT==2):
      parallel_fm3d(velinv)
      copyfile(ttim,mtrav)
      copyfile(ttim,rtrav)
      os.system(resid+'> '+resout)
    if(JOINT==1 or (JOINT==2 and velinv ==1)):
      os.system(surfdirect+' > '+diagnossurf)
      copyfile(ttimessurf,mtravsurf)
  
    try:
        with open(itn,'r') as inviter:
            ITER = int(inviter.readline())
    except:
        ITER = 1
##
    #for ii in range(1,numfiles+1):
    #  os.system(generatemesh+ '  '+str(ii) +' >> '+gmeshout)
# uncomment this
    #os.system(generatemesh+ '  '+str(nomesh) +' >> '+gmeshout)
  
  # Now begin a loop to iteratively apply subspace inversion
  # and FMM
  #
  #os.system(generatemesh+' > '+gmeshout+' '+nomesh)
  #ITER = 1
  print('iterations:',NI)
  while ITER<=NI:
      if velinv == 1:
     # begin random mesh here
         print('interation: '+str(ITER)+' velocity inversion')
         #os.system(subspaceproj+ ' '+nomesh+' > '+subpout)
         #for inet in range(1,numfiles+1):
         #  os.system(subspaceproj+ ' '+nomesh+' '+str(inet)+' > '+subpout)
         #  os.system(backproj_ind+' '+str(inet)+' '+nomesh+' > '+ backprojout)
         for thidx in range(numfiles/nthread):
           threads = []
           for ii in range(1,nthread+1):
             t = parallelrunsub(thidx*nthread+ii)
             t.start()
             threads.append(t)
           
           for t in threads:
             t.join()
 
         copyfile('receivers.in','receiversref.in')
         os.system(combine_sub+' '+str(numfiles)+' >> '+combineout)

         #copyfile('sourcesnew.in','sources.in')
         #copyfile('stimesnew.dat','stimes.dat')
      else:
         copyfile('receivers.in','receiversref.in')
         print('interation: '+str(ITER)+' relocation')
         os.system(inv)
  
      #runid = sys.argv[1]
      iterdir = runid+'/iteration'+str(ITER).zfill(3)
      if not os.path.exists(iterdir):
      	os.makedirs(iterdir)
      files = ['interfaces.in', 'vgrids.in','vgrids_std.in', 'sources.in',\
  	      'mtimes.dat']
      for line in files:
          copyfile(line,iterdir+'/'+line)

      with open(invinput,'r') as fp:
          data = fp.readlines()
      data[22] = lines[ITER+5].split('|')[0].strip()+'\n'
      data[23] = lines[ITER+5].split('|')[1].strip()+'\n'
      data[24] = lines[ITER+5].split('|')[2].strip()+'\n'
      data[30] = lines[3]
      velinv = int(data[22].split()[0])
      with open(invinput,'w') as fp:
          fp.writelines(data)

# uncomment this for iterations
  
      if(JOINT==0 or JOINT==2):
        os.system(frech)
        parallel_fm3d(velinv)
        copyfile(ttim,mtrav)
        os.system(resid+'>> '+resout)
      if(JOINT==1 or (JOINT==2 and velinv ==1)):
        os.system(surfdirect+' > '+diagnossurf)
        copyfile(ttimessurf, mtravsurf)
      ITER=ITER+1
  
  copyfile('vgrids.in',runid+'/vgrids.in')
  copyfile('residuals.dat',runid+'/residuals.dat')
  copyfile('tomo3d_rm.in',runid+'/tomo3d_rm.in')
  end_time = datetime.now()
  print 'Duration: {}'.format(end_time-start_time)
  copyfile('output_sjfz_syn.txt',runid+'/screenout.txt')
  #move('tempdata','tempdata'+str(nomesh))
  return


if __name__ == "__main__":
  import sys
  runid = sys.argv[1]
  nomesh = sys.argv[2]
  if not os.path.isdir('tempdata'):
    os.makedirs('tempdata')
  print('runid:'+runid+' no of mesh: '+nomesh)
  run_randmesh(runid,nomesh)  
