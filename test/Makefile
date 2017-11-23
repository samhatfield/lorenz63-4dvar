# Default compiler
FC = gfortran

default: all

all: test

# Main target: main executable
test: test.o ../params.o ../lorenz63.o ../assim.o
	$(FC) -o $@ $^

# Dependencies
test.o: ../params.o ../lorenz63.o ../assim.o

# Build rules
%.o: %.f90
	$(FC) -c $< -o $(basename $<).o


.PHONY: clean
clean:
	rm -f *.o *.mod main
