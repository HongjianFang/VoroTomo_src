cmd = invert3d
fc = gfortran
iflag = -I/home/zp/local/include
lflag = -L/home/zp/local/lib -lnetcdff
fflags = -O2 -ffixed-line-length-none -ffloat-store\
	 -W -fbounds-check -m64 -mcmodel=medium
fsrcs = lsmrDataModule.f90 lsmrblasInterface.f90 \
	lsmrblas.f90 lsmrModule.f90 aprod.f90 \
	invert3d.f90
objs = $(fsrcs:%.f90=%.o)
all:$(cmd)
$(cmd):$(objs)
	$(fc) $(fflags) $(lflag) $^ -o $@ 
%.o:%.f90
	$(fc) $(fflags) $(iflag) -c $< -o $@

clean:
	rm *.o *.mod $(cmd)
