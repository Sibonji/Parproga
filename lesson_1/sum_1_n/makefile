NUM ?= 4
C_TARGET ?= sum_1_n.c
O_TARGET ?= sum_1_n

compile:
	mpicc $(C_TARGET) -o O_TARGET

run:
	mpirun -np $(NUM) --use-hwthread-cpus O_TARGET

all: compile run
