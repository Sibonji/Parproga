NUM ?= 4

compile:
	mpicc hello_world.c -o hello_world

run:
	mpirun -np $(NUM) --use-hwthread-cpus hello_world

all: compile run
