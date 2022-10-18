all: search

search: main.o input.o Common.o LSH.o
	g++ -o search main.o input.o Common.o LSH.o

main.o: main.cpp
	g++ -c main.cpp

input.o: input.cpp 
	g++ -c  input.cpp 

Common.o: Common.cpp
	g++ -c Common.cpp 

LSH.o: ../NN/LSH/LSH.cpp
	g++ -c ../NN/LSH/LSH.cpp

clear: 
	rm -r *.o search