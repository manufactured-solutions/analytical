all: main
main: main.cpp nsctpl_fwd.hpp nsctpl.hpp test.py
	test.py > test.hpp
	g++ -O0 -fno-unsafe-math-optimizations -Wall -pedantic main.cpp -o main
clean:
	rm -fv *.o main
