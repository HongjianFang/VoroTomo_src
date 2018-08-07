cmd1 = generatemesh 
cmd2 = subspaceproj
cmd3 = tv_backproj
cmd4 = backproj_ind
cmd5 = genratetruemod 
FC = gfortran
#FC = ifort

iflag = -I/home/hfang/tools/netcdf-gfortran/f90
lflag_sub = -L/home/hfang/tools/netcdf-gfortran/lib -lnetcdff
lflag = -L/home/hfang/tools/lapack-3.8.0 -llapack
fflag1 = -O2 
fflag2 = -O2 -fopenmp 
f90_1 = generatemesh.f90
f90_2 = lsmrDataModule.f90 lsmrblasInterface.f90 \
				lsmrblas.f90 lsmrModule.f90 aprod_new.f90 subspaceproj_omp.f90 
f90_3 =	lsmrDataModule.f90 lsmrblasInterface.f90 \
				lsmrblas.f90 lsmrModule.f90 aprod_new.f90 tv_backproj.f90
f90_4 = aprod_new.f90 backproj_ind.f90
f90_5 = generatetruemod.f90

all:${cmd1} ${cmd2} ${cmd3} ${cmd4} ${cmd5}

objs = $(f90_1:%.f90=%.o)
${cmd1}:${objs}
	${FC} $^ -o $@ ${lflag} ${fflag1}

objs = $(f90_2:%.f90=%.o)
${cmd2}:${objs}
	${FC} $^ -o $@ ${lflag} ${lflag_sub} ${iflag} ${fflag2} 

objs = $(f90_3:%.f90=%.o)
${cmd3}:${objs}
	${FC} $^ -o $@ ${fflag1} 

objs = $(f90_4:%.f90=%.o)
${cmd4}:${objs}
	${FC} $^ -o $@ ${fflag1} 

${cmd5}:${objs}
	${FC} $^ -o $@ ${lflag} ${fflag1}
%.o:%.f90
	${FC} $(fflag2) $(iflag) -c $< -o $@


clean:
	rm ${cmd1} ${cmd2} *.o *.mod

