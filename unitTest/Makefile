all:UT

UT:  unitTest.o Common.o meanCurve.o LSH.o 2dCurves.o 1dCurves.o frechet.o curve.o point.o interval.o simplification.o config.o
	g++ -o UT unitTest.o Common.o meanCurve.o LSH.o 2dCurves.o 1dCurves.o frechet.o curve.o point.o interval.o simplification.o config.o -L/usr/local/lib -Wl,-Bstatic -lcunit -Wl,-Bdynamic

unitTest.o: unitTest.cpp
	g++ -c unitTest.cpp

Common.o: ../Common/Common.cpp
	g++ -c ../Common/Common.cpp

meanCurve.o: ../timeSeries/Clustering/meanCurve.cpp 
	g++ -c ../timeSeries/Clustering/meanCurve.cpp 

LSH.o: ../NN/LSH/LSH.cpp unitTest.cpp
	g++ -c ../NN/LSH/LSH.cpp unitTest.cpp

2dCurves.o: ../timeSeries/Search/2dCurves.cpp
	g++ -c ../timeSeries/Search/2dCurves.cpp

1dCurves.o: ../timeSeries/Search/1dCurves.cpp
	g++ -c ../timeSeries/Search/1dCurves.cpp
frechet.o: ../Fred/src/frechet.cpp
	g++ -c ../Fred/src/frechet.cpp

curve.o: ../Fred/src/curve.cpp 
	g++ -c ../Fred/src/curve.cpp

point.o: ../Fred/src/point.cpp
	g++ -c ../Fred/src/point.cpp

interval.o: ../Fred/src/interval.cpp
	g++ -c ../Fred/src/interval.cpp

simplification.o: ../Fred/src/simplification.cpp
	g++ -c ../Fred/src/simplification.cpp

config.o: ../Fred/src/config.cpp 
	g++ -c ../Fred/src/config.cpp

clean: 
	rm -r *.o search