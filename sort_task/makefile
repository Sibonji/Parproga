NUM ?= 4
C_TARGET ?= bitonic.c
O_TARGET ?= bitonic

compile:
	mpicc $(C_TARGET) -lm -o $(O_TARGET)

run:
	mpirun -np $(NUM) --use-hwthread-cpus $(O_TARGET) 256000000

all: compile run