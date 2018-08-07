CMD = surfdirect
FC = gfortran
iflag = -I/home/hfang/tools/netcdf-gfortran/f90
#iflag = 
lflag = -L/home/hfang/tools/netcdf-gfortran/lib -lnetcdff
#lflag = libnetcdff.a  #-L. -lnetcdff
FFLAGS  = -O3 -ffixed-line-length-none -ffloat-store\
	           -fbounds-check -m64 -mcmodel=medium
F90SRC = mod_surf.f90 delsph.f90 surfdirect.f90\
	         initial_mod.f90 savetonetcdf.f90
FSRC = surfdisp96.f
OBJS = $(F90SRC:%.f90=%.o) $(FSRC:%.f=%.o) CalSurfGnew.o
$(CMD):$(OBJS)
	$(FC) -fopenmp $^ $(lflag) -o $@ 
CalSurfGnew.o: CalSurfGnew.f90
	$(FC) -fopenmp $(FFLAGS) -c $< -o $@
%.o: %.f90
	$(FC) $(FFLAGS) $(iflag) -c $(@F:.o=.f90) -o $@
%.o: %.f
	$(FC) $(FFLAGS) -c $(@F:.o=.f) -o $@
clean:
	rm *.o *.mod $(CMD)
