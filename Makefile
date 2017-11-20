# Default compiler
FC = gfortran

# Main target: main executable
main: main.o params.o lorenz63.o utils.o io.o assim.o
	$(FC) -o $@ $^

# Dependencies
main.o: params.o lorenz63.o utils.o io.o assim.o
lorenz63.o: params.o
utils.o: params.o
io.o: params.o
assim.o: lorenz63.o params.o

# Build rules
%.o: %.f90
	$(FC) -c $< -o $(basename $<).o

.PHONY: clean
clean:
	rm -f *.o *.mod main
