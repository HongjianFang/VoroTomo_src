fm3d_objects = 3dfmlib.o propagate.o rays.o frechet.o visual.o 3dfm_main.o matchref.o teleseismic.o
fm3d_modules = mod_3dfm.o
tt_objects = libsun.o libtau.o ellip.o sphdist.o
nn_objects = nn_subsf.o stack.o
svd_objects = svdlib.o

f90comp = ifort
#fflags = -O2  -traceback -heap-arrays -check all -fp-stack-check
fflags = -O2 -I./aktimes/ -I/usr/include

%.o : %.f90
	$(f90comp) $(fflags) -c  $<

%.o : %.f
	$(f90comp) $(fflags) -c  $<

fm3d : $(fm3d_objects) $(fm3d_modules) $(svd_objects) $(nn_objects) $(tt_objects)
	$(f90comp) $(fflags) -o fm3d  $(fm3d_modules) $(fm3d_objects) $(svd_objects) $(nn_objects) $(tt_objects)


$(fm3d_objects) : $(fm3d_modules)


clean:
	rm $(fm3d_objects) $(fm3d_modules) $(svd_objects) $(nn_objects) $(tt_objects)

