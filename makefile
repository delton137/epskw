
COMP = ifort #gfortran
FFLAGS = -O3 -D_FILE_OFFSET_BITS=64 -D_LARGEFILE_SOURCE
FFTW =  #-L/usr/lib -lfftw3 -lm

default : epskw.x 

debug: FFLAGS += -debug -backtrace -traceback -check bounds
debug: epskw.x
	
epskw.x: m_timer.o four1.o main_stuff.o InputOutput.o truncate_datasets.o chi_k.o correlation_function.o calc_phi.o epskw.o 
	$(COMP) m_timer.o four1.o main_stuff.o InputOutput.o truncate_datasets.o chi_k.o correlation_function.o libxdrf.a calc_phi.o epskw.o  $(FFLAGS) $(FFTW) -o epskw.x 

m_timer.o: m_timer.f90
	$(COMP) -c m_timer.f90 $(FFLAGS) 

four1.o: four1.f
	$(COMP) -c four1.f $(FFLAGS) 

InputOutput.o: InputOutput.f90
	$(COMP) -c InputOutput.f90 $(FFLAGS) 

main_stuff.o: main_stuff.f90
	$(COMP) -c main_stuff.f90 $(FFLAGS) 

truncate_datasets.o: truncate_datasets.f90
	$(COMP) -c truncate_datasets.f90 $(FFLAGS) 

chi_k.o: chi_k.f90
	$(COMP) -c chi_k.f90 $(FFLAGS) 

correlation_function.o: correlation_function.f90
	$(COMP) -c correlation_function.f90 $(FFLAGS) $(FFTW)

calc_phi.o: calc_phi.f90
	$(COMP) -c calc_phi.f90 $(FFLAGS) $(FFTW)

epskw.o: epskw.f90
	$(COMP) -c epskw.f90 $(FFLAGS) 

clean:
	rm -rf *.o
	rm -rf *.mod
