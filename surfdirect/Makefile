CMD = surfdirect
FC = gfortran
FFLAGS  = -O3 -ffixed-line-length-none -ffloat-store\
	          -W -fbounds-check -m64 -mcmodel=medium
F90SRC = CalSurfGnew.f90 mod_surf.f90 delsph.f90 surfdirect.f90\
	         initial_mod.f90
FSRC = surfdisp96.f
OBJS = $(F90SRC:%.f90=%.o) $(FSRC:%.f=%.o)
$(CMD):$(OBJS)
	$(FC) -fopenmp $^ -o $@
%.o: %.f90
	$(FC) $(FFLAGS) -c $(@F:.o=.f90) -o $@
%.o: %.f
	$(FC) $(FFLAGS) -c $(@F:.o=.f) -o $@
clean:
	rm *.o *.mod $(CMD)
