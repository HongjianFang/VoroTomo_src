import numpy as np
from scipy import interpolate
#depo = np.array([-1.5, 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0,10.0,11.0,13.0,16.0,20.0,50.0])
depo = np.arange(0,31,1)
velo = np.array([5.0, 5.5 , 5.566 , 5.631 , 5.697 , 5.763 , 5.829 , 5.894 , 5.96 , 6.026 , 6.091 , 6.157 , 6.223 , 6.289 , 6.354 , 6.42 , 6.486 , 6.551 , 6.617 , 6.683 , 6.749 , 6.814 , 6.88 , 6.946 , 7.011 , 7.077 , 7.143 , 7.209 , 7.274 , 7.34 , 7.406 , 7.471 , 7.537 , 7.8])
dep = np.zeros(len(velo),)
dep[0] = -5
dep[1] = -3
dep[-1] = 50
dep [2:-1] = depo

nrad = 15
top = -3
but = 30
#print dep,velo
depn = np.linspace(top,but,nrad)
f = interpolate.interp1d(dep,velo)
print dep[0],velo[0],velo[0]/1.75
for i in range(len(depn)):
    print depn[i],f(depn[i]),f(depn[i])/1.75
print dep[-1],velo[-1],velo[-1]/1.75
