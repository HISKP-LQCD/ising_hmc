MYFLAGS=-Wall -pedantic -Wno-unused-result -O2
MYLIBS=-lm -lgsl -lgslcblas -lfftw3
MYFILES=hmc init lattice aux evolve meas

all: ising_hmc

# Take file identifiers and change empty prefix by ising_ % empty suffix by .o
# Compile all dependencies $^ to target $@
ising_hmc: $(MYFILES:%=ising_%.o)
	gcc $(MYFLAGS) -o $@ $^ $(MYLIBS)

# For every object .o check for .c and .h file with same name
# Compile first dependency $< (the .c file)
%.o: %.c %.h Makefile
	gcc $(MYFLAGS) -c $<

clean:
	rm -f *.o ising_hmc

distclean:
	rm -f *.o *.exe *~ tags ising_hmc
