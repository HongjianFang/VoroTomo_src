def randomchoosing():
#if 1:
  import random
  import numpy as np
  
  with open('abs.dat') as f:
  	ndata = int(f.readline())
  	eqidxsrc = np.zeros(ndata,dtype=int)
  	eqidxrc = np.zeros(ndata,dtype=int)
  	eqidxabst = np.zeros(ndata,)
  	for ii in xrange(ndata):
  		ln = f.readline().split()
  		eqidxsrc[ii] = int(ln[0])-1
  		eqidxrc[ii] = int(ln[1])-1
  		eqidxabst[ii] = float(ln[2])
  		#eqidx[ii,2] = int(ln[2])
  
  nqall = int(np.max(eqidxsrc)+1)
  nnets = 100
  
  eveloc = np.loadtxt('eveloc.dat')
  lon_w = eveloc[:,1].min()
  lon_e = eveloc[:,1].max()
  lat_s = eveloc[:,0].min()
  lat_n = eveloc[:,0].max()
  dep_u = eveloc[:,2].min()
  dep_d = eveloc[:,2].max()
  ncell = 800
  #print lon_w,lon_e,lat_s,lat_n,dep_u,dep_d
  for inet in range(1,nnets+1):
    #print inet	  
    rlat = np.random.rand(ncell,)
    rlat = lat_s+(lat_n-lat_s)*rlat
    rlat = np.pi/2-np.deg2rad(rlat)
    rlon = np.random.rand(ncell,)
    rlon = lon_w+(lon_e-lon_w)*rlon
    rlon = np.deg2rad(rlon)
    rdep = np.random.rand(ncell,)
    rdep = dep_u+(dep_d-dep_u)*rdep
    rdep = 6371-rdep
    
    xpts = rdep*np.sin(rlat)*np.cos(rlon)
    ypts = rdep*np.sin(rlat)*np.sin(rlon)
    zpts = rdep*np.cos(rlat)
    
    nev = len(eveloc)
    eveidx = np.zeros(nev,dtype=int)
    for ii in range(nev): 
  	evlat = np.pi/2-np.deg2rad(eveloc[ii,0])
  	evlon = np.deg2rad(eveloc[ii,1])
  	evdep = 6371-eveloc[ii,2]
  	x_eve = evdep*np.sin(evlat)*np.cos(evlon)
  	y_eve = evdep*np.sin(evlat)*np.sin(evlon)
  	z_eve = evdep*np.cos(evlat)
  	dis = (xpts-x_eve)**2+(ypts-y_eve)**2+(zpts-z_eve)**2
  	eveidx[ii] = np.argmin(dis)
  
    seveidx = np.unique(eveidx)
    #print len(seveidx)
    nchoose = int(len(seveidx)*0.8) 
    choose = np.zeros(nchoose,dtype='int')
  
  #if 1:
    #choose = np.zeros(nchoose,dtype=int) 
    randeveidx = random.sample(xrange(len(seveidx)),nchoose)
    #print seveidx[randeveidx]
    for iev in range(nchoose):
  	  nev_cell = np.where(eveidx==seveidx[randeveidx[iev]])[0]
  	  choose[iev] = random.sample(nev_cell,1)[0]
  
    #choose = random.sample(xrange(nqall),nchoose)
    choosenew = []
    eqidxnew = []
    phasenew = []
    abstnew = []
    for jj in range(nchoose):
      tmp = np.where(eqidxsrc == choose[jj])[0]
      choosenew = np.hstack([choosenew,tmp]) 
      #print len(tmp)
    finalidx = np.array(np.sort(choosenew),int)
    choosenew = finalidx+1
    eqidxsrcnew = eqidxsrc[finalidx]
    eqidxrcnew = eqidxrc[finalidx]
    phasenew = np.ones(len(choosenew),)#eqidx[finalidx,1]
    abstnew = eqidxabst[finalidx]

    with open('otimes.dat') as f:
	    ndata = int(f.readline())
	    srcallidx = np.zeros(ndata,dtype=int)
	    rcallidx = np.zeros(ndata,dtype=int)
	    for ii in range(ndata):
		    ln = f.readline().split()
		    srcallidx[ii] = int(ln[6])
		    rcallidx[ii] = int(ln[7])
  
    with open('tempdata/randomchoose'+str(inet)+'.dat','w') as f:
	    f.write(str(len(choosenew))+'\n')
   	    for ii in range(len(choosenew)):
	        idxsrc = np.where(rcallidx==eqidxrcnew[ii])[0]
	        idxrc = np.where(srcallidx[idxsrc]==eqidxsrcnew[ii])[0]
	        f.write('%d %d %7.4f \n' % (idxsrc[idxrc]+1,phasenew[ii],abstnew[ii]))
    #print 'finishing random picking events'
    #choosenew = np.hstack([len(choosenew),choosenew])
    #phasenew = np.hstack([len(choosenew),phasenew])
    #abstnew = np.hstack([len(abstnew),abstnew])
    #tmp,eqidxnew_idx = np.unique(np.array(eqidxnew,dtype=int),return_inverse=True)
    #eqidxnew = np.sort(np.unique(np.array(eqidxnew,dtype=int)))
    #print ii,len(np.unique(eqidxnew)),len(np.unique(eqidxnew_idx))
    #eqidxnew_idx = np.hstack([len(choosenew),eqidxnew_idx])
    #choosenew = np.vstack([choosenew,eqidxnew_idx,phasenew,abstnew])
#    with open('tempdata/randomchoose'+str(inet)+'.dat','w') as f:
#	    f.write(str(len(choosenew))+'\n')
#	    for ii in range(len(choosenew)):
#		    f.write('%d %d %d %7.4f \n' % (choosenew[ii],eqidxnew_idx[ii],phasenew[ii],abstnew[ii]))
    #np.savetxt('tempdata/randomchoose'+str(inet)+'.dat',choosenew.T)
    #eqidxnew = np.hstack([len(eqidxnew),eqidxnew])
    #np.savetxt('tempdata/updatesrc'+str(inet)+'.dat',eqidxnew,fmt='%d')
