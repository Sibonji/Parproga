C_TARGET   ?= $(O_TARGET).c
O_TARGET   ?= merge

THREAD_NUM ?= 4
ARRAY_SIZE ?= 256

default: compile run

compile:
	mpicc $(C_TARGET) -lm -o $(O_TARGET)

run:
	mpirun -np $(THREAD_NUM) --use-hwthread-cpus $(O_TARGET) $(ARRAY_SIZE)