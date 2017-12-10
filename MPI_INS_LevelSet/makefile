

CC = mpicc
CPPFLAGS = -fopenmp
OPT = -o3

SRC = main.c
OBJ = main.o

all: main
	@echo "-------------------------"
	@echo "Compilation Successful"
	@echo "-------------------------"

main: main.c
	$(CC) $(CPPFLAGS) -c main.c -lm
	$(CC) $(CPPFLAGS) main.o -o main -lm

clean:
	rm -rf *.o main output
