#!/usr/bin/env python
########################################################################################
# Script to run joint inversion with body wave and surface wave data
# Hongjian Fang @ USTC 2016.11.11
########################################################################################
#
import os
from shutil import copyfile
import sys
#
# forward
fmm='fm3d'
frech='frechgen'
ttim='arrivals.dat'
cvg='vgrids.in'
cig='interfaces.in'
csl='sources.in'
diagnos='fm3dlog.out'
# surface wave
#########################################################
surfdirect='surfdirect'
mtravsurf='mtimessurf.dat'
ttimessurf='ttimessurf.dat'
diagnossurf='surfdirect.out'
#
# inversion 
inv='invert3d'
invinput = 'invert3d.in'
itn='inviter.in'
mtrav='mtimes.dat'
rtrav='rtimes.dat'
ivg='vgridsref.in'
iig='interfacesref.in'
isl='sourcesref.in'
#
# residuals
resid='residuals'
resout='residuals.dat'
#
###############################
ifile='tomo3d.in'
with open(ifile, 'r') as fp:
    lines = fp.readlines()
NI = int(lines[0].strip())
BGITER = int(lines[1].strip())
BINV = int(lines[2].strip())
JOINT = int(lines[3].strip())
velinv = 1
###############################
###############################
#
# If necessary, copy the initial velocity,
# interface and source files to the current files
# read by fm3d. Generate the frechet.in file and
# set the iteration number to 1.
#
from datetime import datetime
start_time = datetime.now()
if BGITER == 0:
    copyfile(ivg,cvg)
    copyfile(iig,cig)
    copyfile(isl,csl)
    ITER=1
    with open(itn, 'w') as fp:
        fp.write(str(ITER))
    with open(invinput,'r') as fp:
        data = fp.readlines()
    data[22] = lines[ITER+3].split('|')[0].strip()+'\n'
    data[23] = lines[ITER+3].split('|')[1].strip()+'\n'
    data[24] = lines[ITER+3].split('|')[2].strip()+'\n'
    data[30] = lines[3]
    velinv = int(data[22].split()[0])
    with open(invinput,'w') as fp:
        fp.writelines(data)

    os.system(frech)
#
#
# Run FMM
#
if BGITER==0 or (BGITER==1 and BINV==1 ):
  if(JOINT==0 or JOINT==2):
    import para
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
#
# Now begin a loop to iteratively apply subspace inversion
# and FMM
#
while ITER<=NI:
    os.system(inv)

    runid = sys.argv[1]
    iterdir = runid+'/iteration'+str(ITER).zfill(3)
    if not os.path.exists(iterdir):
    	os.makedirs(iterdir)
    files = ['interfaces.in', 'vgrids.in', 'sources.in',\
	      'mtimes.dat','mtimessurf.dat','invert3d.in','inviter.in','vgrids.invpvs']
    for line in files:
        copyfile(line,iterdir+'/'+line)

    with open(invinput,'r') as fp:
        data = fp.readlines()
    data[22] = lines[ITER+4].split('|')[0].strip()+'\n'
    data[23] = lines[ITER+4].split('|')[1].strip()+'\n'
    data[24] = lines[ITER+4].split('|')[2].strip()+'\n'
    data[30] = lines[3]
    velinv = int(data[22].split()[0])
    with open(invinput,'w') as fp:
        fp.writelines(data)

    if(JOINT==0 or JOINT==2):
      os.system(frech)
      if 'para' in sys.modules:
         reload(para)
      else:
	 import para
      copyfile(ttim,mtrav)
      os.system(resid+'>> '+resout)
    if(JOINT==1 or (JOINT==2 and velinv ==1)):
      os.system(surfdirect+' > '+diagnossurf)
      copyfile(ttimessurf, mtravsurf)
    ITER=ITER+1
    COUNT=ITER

    with open(itn, 'w') as fp:
        fp.write(str(COUNT))

copyfile('vgrids.in',runid+'/vgrids.in')
copyfile('vgrids.invpvs',runid+'/vgrids.invpvs')
copyfile('tomo3d.in',runid+'/tomo3d.in')
copyfile('invert3d.in',runid+'/invert3d.in')
copyfile('residuals.dat',runid+'/residuals.dat')
end_time = datetime.now()
print 'Duration: {}'.format(end_time-start_time)
copyfile('screenout.txt',runid+'/screenout.txt')
