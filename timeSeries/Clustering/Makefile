CC =  g++  -std=c++14
objects=main.o Lloyd.o LSH.o Common.o Proj.o KMeansRangeSearch.o meanCurve.o \
		KMeansCurves.o 2dCurves.o KMeansRangeSearchCurves.o 1dCurves.o \
		frechet.o	curve.o	interval.o simplification.o point.o config.o
cflags =  -g3 -Wall 
program = cluster 


$(program):  $(objects)
	$(CC) -o $(program) $(objects)


2dCurves.o:../Search/2dCurves.cpp 
	$(CC) $(cflags) -c ../Search/2dCurves.cpp

KMeansCurves.o:KMeansCurves.cpp 
	$(CC) $(cflags) -c KMeansCurves.cpp
meanCurve.o:meanCurve.cpp 
	$(CC) $(cflags) -c meanCurve.cpp
KMeansRangeSearch.o: KMeansRangeSearch.cpp 
	$(CC) $(cflags) -c KMeansRangeSearch.cpp
KMeansRangeSearchCurves.o:KMeansRangeSearchCurves.cpp 
	$(CC) $(cflags) -c KMeansRangeSearchCurves.cpp
Lloyd.o: Lloyd.cpp
	$(CC) $(cflags) -c Lloyd.cpp
LSH.o: ../../NN/LSH/LSH.cpp
	$(CC) $(cflags) -c ../../NN/LSH/LSH.cpp
Proj.o: ../../NN/Proj/Proj.cpp
	$(CC) $(cflags) -c ../../NN/Proj/Proj.cpp
Common.o: ../../Common/Common.cpp
	$(CC) $(cflags) -c  ../../Common/Common.cpp

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

1dCurves.o:../Search/1dCurves.cpp
	$(CC) $(cflags) -c ../Search/1dCurves.cpp
run: $(program)
	./$(program)

debug: $(program)
	gdb ./$(program)



clean: 
	rm -f $(program)
	rm -f $(objects)