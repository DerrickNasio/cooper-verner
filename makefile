# Usage:
# make        # compile all binary
# make clean  # remove ALL binaries and objects

.PHONY = clean

CC = gcc  # compiler to use

Verner : main.o tests.o Verner.o
	@echo "Linking..."
	${CC} main.o tests.o Verner.o -lc -lm -o Verner
	@echo "Done."

main.o : main.c tests.h Verner.h
	@echo "Compiling main method..."
	${CC} -c -g main.c -o main.o

tests.o : tests.c tests.h
	@echo "Compiling test library..."
	${CC} -c -g tests.c -o tests.o

Verner.o : butcher_tableau.c dense_output.c Verner.c Verner.h
	@echo "Compiling integrator library..."
	${CC} -c -g butcher_tableau.c dense_output.c Verner.c -o Verner.o

clean :
	@echo "Cleaning up..."
	\rm main.o tests.o Verner.o Verner
	@echo "Done."
