#!/bin/bash
#SBATCH -N 1
#SBATCH -n 5 
#SBATCH -p sched_mit_hill
#SBATCH --time=5:00:00
#SBATCH --mem=50000
#SBATCH --exclusive
#SBATCH -o output_mogul.txt

module load gcc/4.8.4
gfortran --version
gcc --version
which gfortran
which gcc
#cd /home/hfang/tools/netcdf/netcdf-4.4.1
#export CC=gcc
#./configure --prefix=/home/hfang/tools/netcdf/netcdf-C --disable-dap --disable-netcdf-4
#make check
#make install
#cd /home/hfang/tools/netcdf/netcdf-fortran-4.4.4
#export CC=gcc
#export FC=gfortran
#export CPPFLAGS=-I/home/hfang/tools/netcdf/netcdf-C/include
#export LDFLAGS=-L/home/hfang/tools/netcdf/netcdf-C/lib
#./configure --prefix=/home/hfang/tools/netcdf/netcdf-Fortran
#make check
#make install
cd /pool001/hfang/vorotomo/VoroTomo_src/fm3d_src
make -f Makefile.eofe7 clean
make -f Makefile.eofe7
cd ../mogultomo
#cd ../data_Malcolm_absrtt

#LD_LIBRARY_PATH = /home/hfang/tools/netcdf/netcdf-Fortran/lib:$LD_LIBRARY_PATH
#export LD_LIBRARY_PATH
echo $LD_LIBRARY_PATH
#export OMP_NUM_THREADS=10
#python VoroTomo.py randmesh_500_abs_real5 96
#python VoroTomo_rtt.py randmesh_500_rtt_real_absini2 96
python VoroTomo_rtt.py randmesh_500 100
