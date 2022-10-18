CC =  g++  -std=c++14
objects = LSH.o test.o Common.o 2dCurves.o OutputFunct.o Proj.o input.o 1dCurves.o \
	frechet.o	curve.o	interval.o simplification.o point.o config.o


cflags =  -g3 
program = search

$(program):  $(objects)
	$(CC) -o $(program) $(objects)

1dCurves.o:1dCurves.cpp
	$(CC) $(cflags) -c 1dCurves.cpp
OutputFunct.o:OutputFunct.cpp
	$(CC) $(cflags) -c OutputFunct.cpp

2dCurves.o:2dCurves.cpp 
	$(CC) $(cflags) -c 2dCurves.cpp
test.o: test.cpp
	$(CC) $(cflags) -c test.cpp

LSH.o: ../../NN/LSH/LSH.cpp
	$(CC) $(cflags) -c ../../NN/LSH/LSH.cpp
	$(CC) $(cflags) -c ../../NN/Proj/Proj.cpp
Common.o: ../../Common/Common.cpp
	$(CC) $(cflags) -c  ../../Common/Common.cpp

input.o: ../../Common/input.cpp
	$(CC) $(cflags) -c ../../Common/input.cpp
frechet.o: ../../Fred/src/frechet.cpp
	g++ -c ../../Fred/src/frechet.cpp

curve.o: ../../Fred/src/curve.cpp 
	g++ -c ../../Fred/src/curve.cpp

point.o: ../../Fred/src/point.cpp
	g++ -c ../../Fred/src/point.cpp

interval.o: ../../Fred/src/interval.cpp
	g++ -c ../../Fred/src/interval.cpp

simplification.o: ../../Fred/src/simplification.cpp
	g++ -c ../../Fred/src/simplification.cpp

config.o: ../../Fred/src/config.cpp 
	g++ -c ../../Fred/src/config.cpp





run: $(program)
	./$(program)

debug: $(program)
	gdb ./$(program)

reset:
	rm -f test.o

clean: 
	rm -f $(program)
	rm -f $(objects)
