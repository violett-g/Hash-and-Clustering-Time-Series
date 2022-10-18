CC =  g++  -std=c++11
objects = Proj.o test.o Common.o LSH.o

cflags =  -g3 -Wall
program = cube

$(program):  $(objects)
	$(CC) -o $(program) $(objects)



test.o: test.cpp
	$(CC) $(cflags) -c test.cpp

Proj.o: Proj.cpp
	$(CC) $(cflags) -c Proj.cpp
LSH.o: ../LSH/LSH.cpp
	$(CC) $(cflags) -c ../LSH/LSH.cpp
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