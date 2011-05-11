all: main

test.hpp: test.py soln.py forcing.py
	test.py | tee test.hpp

main: main.cpp nsctpl_fwd.hpp nsctpl.hpp test.hpp
	g++ -O0 -fno-unsafe-math-optimizations -Wall -pedantic main.cpp -o main

clean:
	rm -fv *.o main test.hpp soln.c soln.h
