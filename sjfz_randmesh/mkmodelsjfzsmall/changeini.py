import numpy as np
with open('vgridsref.in','r') as vref:
   ln = vref.readline().split()
   vtp = int(ln[1])
   dim = vref.readline().split()
   nr = int(dim[0])
   nlat = int(dim[1])
   nlon = int(dim[2])
   line3 = vref.readline()
   line4 = vref.readline()
   vp = np.zeros((nlon,nlat,nr)) 
   vs = np.zeros((nlon,nlat,nr)) 
   for k in range(nr):
    for j in range(nlat):
      for i in range(nlon):
        vp[i,j,k] = float(vref.readline().split()[0])
   for k in range(nr//3,nr*2//3):
    for j in range(nlat//3,nlat*2//3):
      for i in range(nlon//3,nlon*2//3):
        vp[i,j,k] = vp[i,j,k]-0.4

   vref.readline()
   vref.readline()
   vref.readline()
   for k in range(nr):
    for j in range(nlat):
      for i in range(nlon):
        vs[i,j,k] = float(vref.readline().split()[0])
   for k in range(2,nr//2):
    for j in range(nlat//3,nlat*2//3):
      for i in range(nlon//3,nlon*2//3):
        vs[i,j,k] = vp[i,j,k]/1.6
   for k in range(nr//2+1,nr):
    for j in range(nlat//3,nlat*2//3):
      for i in range(nlon//3,nlon*2//3):
        vs[i,j,k] = vp[i,j,k]/1.9

with open('vgridsnew.in','w') as vn:
 vn.write(' '.join(['1','2','\n']))
 vn.write(' '.join([str(nr),str(nlat),str(nlon),'\n']))
 vn.write(line3)
 vn.write(line4)
 for k in range(nr):
    for j in range(nlat):
      for i in range(nlon):
 	vn.write(' '.join([str(vp[i,j,k]),'\n']))
 vn.write(' '.join([str(nr),str(nlat),str(nlon),'\n']))
 vn.write(line3)
 vn.write(line4)
 for k in range(nr):
    for j in range(nlat):
      for i in range(nlon):
 	vn.write(' '.join([str(vs[i,j,k]),'\n']))
