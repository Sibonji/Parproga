NUM ?= 1
C_TARGET_CL ?= client.c
O_TARGET_CL ?= client

C_TARGET_SE ?= server.c
O_TARGET_SE ?= server
PORT ?= 

compile_se:
	mpicc $(C_TARGET_SE) -o $(O_TARGET_SE)

run_se:
	mpiexec -np $(NUM) --use-hwthread-cpus $(O_TARGET_SE) $(PORT)

compile_cl:
	mpicc $(C_TARGET_CL) -o $(O_TARGET_CL)

run_cl:
	mpiexec -np $(NUM) --use-hwthread-cpus $(O_TARGET_CL) $(PORT)

all_cl: compile_cl run_cl

all_se: compile_se run_se