# Makefile

FC      = gfortran
FFLAGS  = -O3 -mcmodel=medium -llapack -lblas
LIBS    = -llapack -lblas
INCLUDE =

OBJS    = m_gridsse.o m_dc3d_driver.o dc3d.o m_rtrend.o m_sort.o m_inv.o

.SUFFIXES:
.SUFFIXES:.f90 .o .f


gridsse: $(OBJS) gridsse.f90
	$(FC) $(FFLAGS) $^ -o $@ $(LIBS) $(INCLUDE)


.f90.o:
	$(FC) $(FFLAGS) -c $< -o $@ $(LIBS) $(INCLUDE)
m_sort.o        : $(@:.o=.f90)
m_rtrend.o      : $(@:.o=.f90)
m_inv.o         : $(@:.o=.f90)
dc3d.o          : $(@:.o=.f)
m_dc3d_driver.o : $(@:.o=.f90) dc3d.o
m_gridsse.o     : $(@:.o=.f90) dc3d.o m_dc3d_driver.o m_rtrend.o m_sort.o m_inv.o

clean:
	rm -f *~ *.o *.mod

