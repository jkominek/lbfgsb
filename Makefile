FC = gfortran

FFLAGS = -O -Wall -fbounds-check -g -Wno-uninitialized 

LBFGSB  = lbfgsb.f
LINPACK = linpack.f
BLAS    = blas.f
TIMER   = timer.f

all :  lbfgsb

lbfgsb :  $(LBFGSB) $(LINPACK) $(BLAS) $(TIMER)
	$(FC) $(FFLAGS) $(LBFGSB) $(LINPACK) $(BLAS) $(TIMER) -fPIC -shared -o lbfgsb.so
