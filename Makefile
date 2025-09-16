CC=g++
CFLAGS= -O0 -mfpmath=sse -fstack-protector-all -g -W -Wall -Wextra -Wunused -Wcast-align -Werror -pedantic -pedantic-errors -Wfloat-equal -Wpointer-arith -Wformat-security -Wmissing-format-attribute -Wformat=1 -Wwrite-strings -Wcast-align -Wno-long-long -Woverloaded-virtual -Wnon-virtual-dtor -Wcast-qual -ffast-math -Wno-suggest-attribute=format

a.out: main.o initialize_matrix.o progonka_solver.o log.o
	$(CC) -g main.o initialize_matrix.o progonka_solver.o log.o -o a.out

main.o : main.cpp
	$(CC) $(CFLAGS) -c main.cpp -o main.o 

initialize_matrix.o: initialize_matrix.cpp
	$(CC) $(CFLAGS) -c initialize_matrix.cpp -o initialize_matrix.o 

progonka_solver.o : progonka_solver.cpp
	$(CC) $(CFLAGS) -c progonka_solver.cpp -o progonka_solver.o

log.o : log.cpp
	$(CC) $(CFLAGS) -c log.cpp -o log.o

clean:
	rm -f *.o a.out