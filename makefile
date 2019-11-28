all: test

test: test.o
	g++ test.o -o test

test.o: test.cpp
	g++ -c test.cpp

clear:
	rm *.o *.out test out.txt

plot:
	gnuplot figure.plt