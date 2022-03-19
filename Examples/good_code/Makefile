CC = g++
CFLAGS = -Wall -O3 --std=c++11 `root-config --cflags`
SOURCES=main.cpp function.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=main


main : function.o main.o
	$(CC) function.o main.o -o main $(CFLAGS)

function.o : function.cpp function.h
	$(CC) -c function.cpp -o function.o $(CFLAGS)

main.o : main.cpp function.h
	$(CC) -c main.cpp -o main.o $(CFLAGS)

clean :
	rm *.o main
