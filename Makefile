EXECUTABLE=xsect.exe
SOURCES=types_module.f90 binomial_module.f90 \
abrasion.f90 density.f90 lgwt.f90 main.f90
OBJECTS=$(SOURCES:.f90=.o)
LDFLAGS=

FC=gfortran
CFLAGS=-g -O3 -ffree-line-length-none -Wall -Wextra -Wconversion \
-pedantic -ffpe-trap=zero,overflow,invalid \
-finit-local-zero -fbounds-check -fcheck-array-temporaries

all: $(SOURCES) $(EXECUTABLE)
	
$(EXECUTABLE): $(OBJECTS) 
	$(FC) $(LDFLAGS) $(OBJECTS) -o $@

%.o: %.f90
	$(FC) -c $(CFLAGS) $< -o $@

clean:
	rm -f *.o *.mod $(EXECUTABLE)
