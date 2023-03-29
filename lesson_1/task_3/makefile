NUM ?= 4
C_TARGET ?= task_3.c
O_TARGET ?= task_3

compile:
	mpicc $(C_TARGET) -o O_TARGET

run:
	mpirun -np $(NUM) --use-hwthread-cpus O_TARGET

all: compile run
