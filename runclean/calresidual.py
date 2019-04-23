#/* -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.
#
# File Name : residual.py
#
# Purpose :
#
# Creation Date : 16-01-2019
#
# Last Modified : Tue 19 Feb 2019 08:56:37 AM EST
#
# Created By : Hongjian Fang: hfang@mit.edu
#
#_._._._._._._._._._._._._._._._._._._._._.*/
def calresidual():
#if 1:
   import pandas
   import random
   import numpy as np
   #pd = pandas.read_csv('differential_P.txt',sep=' ',header=None)
   pd = pandas.read_table('differential_P.txt',delim_whitespace=True,header=None)
   pdsyn=pandas.read_table('arrivals.dat',delim_whitespace=True,header=None,names=('i1','i2','i3','i4','i5','i6','i7'))
   pdobs=pandas.read_table('otimes.dat',delim_whitespace=True,header=0,names=('i1','i2','i3','i4','i5','i6','i7','i8'))
   
   nd = len(pd)
   nchoose = 10000
   idx = random.sample(xrange(nd),nchoose)
   stidx = pdobs['i8']
   obs = np.zeros(nchoose,)
   syn = np.zeros(nchoose,)
   for jj in range(nchoose):
   	ii = idx[jj]
   	obs[jj] = pd.loc[ii,3]
   	evid1 = pd.loc[ii,0]-1
   	evid2 = pd.loc[ii,1]-1
   	staid = pd.loc[ii,2]-1
   	idxst = np.where(stidx==staid)[0]
   	evidx = pdobs['i7'][idxst]
   	idx1 = np.where(evidx==evid1)[0]
   	idx2 = np.where(evidx==evid2)[0]
   	if len(idx1)>0 and len(idx2)>0:
   	  evetime1 = float(pdsyn['i5'][idxst[idx1]])
   	  evetime2 = float(pdsyn['i5'][idxst[idx2]])
   	  syn[jj] = evetime2-evetime1
   err = obs-syn
   print 'mean and std of rtt:',err.mean(),err.std()
   weight = 1.0/(1+0.05*(err**2*5.0))
   errweight = err*weight
   print 'mean and std of rtt after weighting:',errweight.mean(),errweight.std()
