cmd1 = generatemesh 
cmd2 = subspaceproj
cmd3 = subspaceproj_rtt
cmd23 = subspaceproj_absrtt
#cmd3 = tv_backproj
cmd4 = backproj_ind
#cmd5 = genratetruemod 
cmd6 = combine_sub
cmd7 = relocation
cmd8 = relocation_rtt
cmd78 = relocation_absrtt
FC = gfortran
#FC = ifort

#iflag = -I/home/seis/ping/Documents/netcdf-fortran-4.4.4/include
#iflag = -I/home/hfang/tools/netcdf-gfortran/f90
iflag = -I/home/hfang/tools/netcdf/netcdf-Fortran/include
lflag_sub = -L/home/hfang/tools/netcdf/netcdf-Fortran/lib -lnetcdff
#lflag_sub = -L/home/hfang/tools/netcdf-gfortran/lib -lnetcdff
lflag = -L/home/hfang/tools/lapack-3.8.0 -llapack
#lflag_sub = -L/home/seis/ping/Documents/netcdf-fortran-4.4.4/lib -lnetcdff
#lflag = -L/etc/alternatives -llapack
fflag1 = -O2 
#fflag2 = -O2 -f#openmp 
f90_1 = generatemesh.f90
f90_2 = lsmrDataModule.f90 lsmrblasInterface.f90 \
	lsmrblas.f90 lsmrModule.f90 aprod_new.f90 subspaceproj.f90 

f90_3 = lsmrDataModule.f90 lsmrblasInterface.f90 \
	lsmrblas.f90 lsmrModule.f90 aprod_new.f90 subspaceproj_rtt.f90 

f90_23 = lsmrDataModule.f90 lsmrblasInterface.f90 \
	lsmrblas.f90 lsmrModule.f90 aprod_new.f90 subspaceproj_absrtt.f90 

f90_4 = aprod_new.f90 backproj_ind.f90
f90_6 = combine_sub.f90
f90_7 = lsmrDataModule.f90 lsmrblasInterface.f90 \
	lsmrblas.f90 lsmrModule.f90 aprod_new.f90 relocation.f90 
f90_8 = lsmrDataModule.f90 lsmrblasInterface.f90 \
	lsmrblas.f90 lsmrModule.f90 aprod_new.f90 relocation_rtt.f90 
f90_78 = lsmrDataModule.f90 lsmrblasInterface.f90 \
	lsmrblas.f90 lsmrModule.f90 aprod_new.f90 relocation_absrtt.f90 

all:${cmd1} ${cmd2} ${cmd3} ${cmd23} ${cmd4} ${cmd6} ${cmd7} ${cmd8} ${cmd78}

objs = $(f90_1:%.f90=%.o)
${cmd1}:${objs}
	${FC} $^ -o $@ ${lflag} ${fflag1}

objs = $(f90_2:%.f90=%.o)
${cmd2}:${objs}
	${FC} $^ -o $@ ${lflag} ${lflag_sub} ${fflag2} 

objs = $(f90_3:%.f90=%.o)
${cmd3}:${objs}
	${FC} $^ -o $@ ${lflag} ${lflag_sub} ${fflag2} 

objs = $(f90_23:%.f90=%.o)
${cmd23}:${objs}
	${FC} $^ -o $@ ${lflag} ${lflag_sub} ${fflag2} 

objs = $(f90_4:%.f90=%.o)
${cmd4}:${objs}
	${FC} $^ -o $@ ${fflag1} 

objs = $(f90_5:%.f90=%.o)
${cmd5}:${objs}
	${FC} $^ -o $@ ${lflag} ${fflag1}

objs = $(f90_6:%.f90=%.o)
${cmd6}:${objs}
	${FC} $^ -o $@ ${lflag} ${fflag1}

objs = $(f90_7:%.f90=%.o)
${cmd7}:${objs}
	${FC} $^ -o $@ ${lflag} ${lflag_sub} ${iflag} ${fflag2} 

objs = $(f90_8:%.f90=%.o)
${cmd8}:${objs}
	${FC} $^ -o $@ ${lflag} ${lflag_sub} ${iflag} ${fflag2} 

objs = $(f90_78:%.f90=%.o)
${cmd78}:${objs}
	${FC} $^ -o $@ ${lflag} ${lflag_sub} ${iflag} ${fflag2} 

%.o:%.f90
	${FC} $(fflag2) $(iflag) -c $< -o $@


clean:
	rm ${cmd1} ${cmd2} *.o *.mod

