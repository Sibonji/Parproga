LAB        := lab
DOPUSK     := dopusk
THREAD_NUM ?= 4
K          ?= 1000
M          ?= 1000

CompileDopusk:
	mpicc $(DOPUSK).c -o $(DOPUSK)

RunDopusk:
	mpiexec -n 2 ./$(DOPUSK) 10000000

CompileLab:
	mpicc $(LAB).c -lm -o $(LAB)

RunLab:
	mpiexec -n $(THREAD_NUM) ./$(LAB) $(K) $(M)

Dopusk: CompileDopusk RunDopusk

Lab: CompileLab RunLab