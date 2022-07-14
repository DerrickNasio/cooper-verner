# Usage:
# make        # compile all binary
# make clean  # remove ALL binaries and objects

.PHONY = clean

CC = gcc  # compiler to use

Verner : main.o tests.o Verner.o
	@echo "Linking object files..."
	${CC} main.o tests.o Verner.o -o Verner
	@echo "Done."

main.o : main.c tests.o Verner.o
	@echo "Compiling main method..."
	${CC} -c main.c

tests.o : tests.c tests.h
	@echo "Compiling test library..."
	${CC} -c tests.c

Verner.o : butcher_tableau.c dense_output.c Verner.c Verner.h
	@echo "Compiling integrator library..."
	${CC} -c butcher_tableau.c dense_output.c Verner.c

clean :
	@echo "Cleaning up..."
	\rm main.o tests.o Verner.o Verner
	@echo "Done."
