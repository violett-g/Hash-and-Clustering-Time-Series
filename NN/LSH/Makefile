CC =  g++  -std=c++11
objects = LSH.o test.o Common.o

cflags =  -g3 -Wall
program = lsh

$(program):  $(objects)
	$(CC) -o $(program) $(objects)



test.o: test.cpp
	$(CC) $(cflags) -c test.cpp

LSH.o: LSH.cpp
	$(CC) $(cflags) -c LSH.cpp
Common.o: ../../Common/Common.cpp
	$(CC) $(cflags) -c  ../../Common/Common.cpp

run: $(program)
	./$(program)

debug: $(program)
	gdb ./$(program)

reset:
	rm -f test.o

clean: 
	rm -f $(program)
	rm -f $(objects)