FC=gfortran
FFLAGS=-fopenmp -fno-unsafe-math-optimizations -frounding-math -fsignaling-nans -O3 -march=native -ffast-math
TARGET:meteotsun

.SUFFIXES:
.SUFFIXES: .f90 .o

mod_params.o: mod_params.f90
mod_grid.o: mod_grid.f90 mod_params.o
mod_bathy.o: mod_bathy.f90 mod_params.o
mod_gauge.o: mod_gauge.f90 mod_params.o
mod_boundary.o: mod_boundary.f90 mod_params.o
mod_press.o: mod_press.f90 mod_params.o
mod_subroutines.o: mod_subroutines.f90 mod_params.o

.f90.o:
	$(FC) $(FFLAGS) -c $< -o $@

meteotsun: mod_params.o mod_grid.o mod_bathy.o mod_gauge.o mod_boundary.o mod_press.o mod_subroutines.o main_tdac.f90
	$(FC) $(FFLAGS) -o $@ $^

clean:
	/bin/rm -f *.o *.mod *~ meteotsun
