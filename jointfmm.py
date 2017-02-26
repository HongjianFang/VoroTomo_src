#!/usr/bin/env python
########################################################################################
# Script to run joint inversion with body wave and surface wave data
# Hongjian Fang @ USTC 2016.11.11
########################################################################################
#
import os
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
    os.system('cp %s %s' % (ivg, cvg))
    os.system('cp %s %s' % (iig, cig))
    os.system('cp %s %s' % (isl, csl))
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
    #os.system(fmm+' > '+diagnos)
    import para
    os.system('cp %s %s' % (ttim, mtrav))
    os.system('cp %s %s' % (ttim, rtrav))
    os.system(resid+'> '+resout)
  if(JOINT==1 or (JOINT==2 and velinv ==1)):
    os.system(surfdirect+' > '+diagnossurf)
    os.system('cp %s %s' % (ttimessurf, mtravsurf))

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

    for line in ['interfaces.in', 'vgrids.in', 'sources.in']:
        os.system('cp %s %s' % (line, line+'.'+str(ITER).zfill(3)))

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
      #os.system(fmm+' > '+diagnos)
      reload(para)
      os.system('cp %s %s' % (ttim, mtrav))
      os.system(resid+'>> '+resout)
    if(JOINT==1 or (JOINT==2 and velinv ==1)):
      os.system(surfdirect+' > '+diagnossurf)
      os.system('cp %s %s' % (ttimessurf, mtravsurf))
    ITER=ITER+1
    COUNT=ITER

    with open(itn, 'w') as fp:
        fp.write(str(COUNT))

end_time = datetime.now()
print 'Duration: {}'.format(end_time-start_time)
