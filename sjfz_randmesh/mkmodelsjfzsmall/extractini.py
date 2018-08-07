#!/usr/bin/env python
# extract initial model from MOD from tomoDD
# writen by Hongjian Fang @USTC Nov 19, 2016
import numpy as np
from scipy import interpolate

lon_ori1 = 241.83
lon_ori2 = 244.62
lat_ori1 = 32.38
lat_ori2 = 34.54
dep1 = -1.5
dep2 = 30.0
nlon = 100
nlat = 80
ndep = 30
pltfig = False
#pltfig = True


with open('modinv_ng.dat','r') as mod:
    lins = mod.readline()
    nx = int(lins.split()[1])
    ny = int(lins.split()[2])
    nz = int(lins.split()[3])
    lon = [float(x) for x in mod.readline().split()]
    lat = [float(x) for x in mod.readline().split()]
    dep = [float(x) for x in mod.readline().split()]
    vel = np.zeros((nx,ny,nz))
    vels = np.zeros((nx,ny,nz))
    for i in range(nz):
    	for j in range(ny):
	   vel[:,j,i] = [float(x) for x in mod.readline().split()]
    for i in range(nz):
    	for j in range(ny):
	   vels[:,j,i] = [float(x) for x in mod.readline().split()]

lonnew = np.zeros(nlon+2,)
latnew = np.zeros(nlat+2,)
depnew = np.zeros(ndep+2,)
lonnew[1:-1] = np.linspace(lon_ori1,lon_ori2,nlon)
lonnew[0] = lonnew[1] - 0.1
lonnew[-1] = lonnew[-2] + 0.1
latnew[1:-1] = np.linspace(lat_ori1,lat_ori2,nlat)
latnew[0] = latnew[1] - 0.1
latnew[-1] = latnew[-2] + 0.1
depnew[1:-1] = np.linspace(dep1,dep2,ndep)
depnew[0] = depnew[1] - 0.1
depnew[-1] = depnew[-2] + 0.1

velt = np.zeros((nx,ny,ndep+2))
velnew = np.zeros((nlon+2,nlat+2,ndep+2))
velts = np.zeros((nx,ny,ndep+2))
velnews = np.zeros((nlon+2,nlat+2,ndep+2))

for j in range(ny):
  for i in range(nx):
	velt[i,j,:] = np.interp(depnew,dep,vel[i,j,:])

for j in range(ny):
  for i in range(nx):
	velts[i,j,:] = np.interp(depnew,dep,vels[i,j,:])

if pltfig:
  from matplotlib import pyplot as plt
  plt.figure()
  plt.imshow(velt[:,:,5])
  plt.colorbar()

for k in range(ndep+2):
   xx,yy = np.meshgrid(lat,lon)
   velt1 = velt[:,:,k]
   f = interpolate.interp2d(lat,lon,velt1,kind='linear')
   velnew[:,:,k] = f(latnew,lonnew)
for k in range(ndep+2):
   xx,yy = np.meshgrid(lat,lon)
   velt1 = velts[:,:,k]
   f = interpolate.interp2d(lat,lon,velt1,kind='linear')
   velnews[:,:,k] = f(latnew,lonnew)

if pltfig:
  plt.figure()
  plt.imshow(velnew[:,:,5])
  plt.colorbar()
  plt.show()

for k in range(ndep+2):
 for j in range(nlat+2):
   for i in range(nlon+2):
     velp =velnew[i,j,ndep+1-k]
     print '%7.4f' % velp 
for k in range(ndep+2):
 for j in range(nlat+2):
   for i in range(nlon+2):
     vels =velnews[i,j,ndep+1-k]
     print '%7.4f' % vels 

