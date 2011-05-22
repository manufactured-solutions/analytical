MAKEFLAGS := -k

all: writeup.pdf check module

writeup.pdf: writeup.tex soln.py forcing.py nsctpl_fwd.hpp nsctpl.hpp
	md5sum soln.py forcing.py nsctpl_fwd.hpp nsctpl.hpp > CHECKSUMS
	pdflatex writeup
	pdflatex writeup

test.hpp: test.py soln.py forcing.py
	test.py | tee test.hpp

check: check.cpp nsctpl_fwd.hpp nsctpl.hpp test.hpp
	g++ -g -O0 -fno-unsafe-math-optimizations -Wall -pedantic check.cpp -o check

module: _nsctpl.so
_nsctpl.so: nsctpl_fwd.hpp nsctpl.hpp nsctpl.i check
	swig -v -Werror -c++ -python nsctpl.i
	python setup.py build_ext --inplace clean

clean:
	rm -fv *.o check test.hpp soln.c soln.h *.aux *.log *.out writeup.pdf
	rm -fv nsctpl.py nsctpl_wrap.cxx

