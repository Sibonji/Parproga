LAB        := lab
NONPARLAB  := non_par_lab
DOPUSK     := dopusk
THREAD_NUM ?= 2
K          ?= 1000
M          ?= $(K)

CompileDopusk:
	mpicc $(DOPUSK).c -o $(DOPUSK)

RunDopusk:
	mpiexec -n 2 ./$(DOPUSK) 10000000

CompileLab:
	mpicc $(LAB).c -lm -o $(LAB)

RunLab:
	mpiexec -n $(THREAD_NUM) ./$(LAB) $(K) $(M)

CompileNonParLab:
	gcc $(NONPARLAB).c -lm -o $(NONPARLAB)

RunNonParLab:
	./$(NONPARLAB) $(K) $(M)

Dopusk: CompileDopusk RunDopusk

Lab: CompileLab RunLab
NonParLab: CompileNonParLab RunNonParLab