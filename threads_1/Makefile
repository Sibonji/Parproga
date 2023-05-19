SUM=sum
HELLO=hello
CELL=cell

THREAD_NUM ?= 2

run_hello:
	mpiCC $(HELLO).cpp -o $(HELLO)
	mpirun -np $(THREAD_NUM) ./$(HELLO) 3

run_sum:
	mpiCC $(SUM).cpp -o $(SUM)
	mpirun -np $(THREAD_NUM) ./$(SUM) 5 4

run_cell:
	mpiCC $(CELL).cpp -o $(CELL)
	mpirun -np $(THREAD_NUM) ./$(CELL) 2