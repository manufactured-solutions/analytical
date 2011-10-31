MAKEFLAGS := -k

all: writeup.pdf check module

writeup.pdf: writeup.tex soln.py forcing.py nsctpl_rholut_fwd.hpp nsctpl_rholut.hpp README
	md5sum soln.py forcing.py nsctpl_rholut_fwd.hpp nsctpl_rholut.hpp > CHECKSUMS
	pdflatex writeup
	pdflatex writeup

test.hpp: test.py soln.py forcing.py
	test.py | tee test.hpp

check: check.cpp nsctpl_rholut_fwd.hpp nsctpl_rholut.hpp test.hpp
	g++ -g -O0 -fno-unsafe-math-optimizations -Wall -pedantic check.cpp -o check

module: _nsctpl_rholut.so
_nsctpl_rholut.so: nsctpl_rholut_fwd.hpp nsctpl_rholut.hpp nsctpl_rholut.i check
	swig -v -Werror -c++ -python nsctpl_rholut.i
	python setup.py build_ext --inplace clean

clean:
	rm -fv *.o check test.hpp soln.c soln.h *.aux *.log *.out writeup.pdf
	rm -fv nsctpl_rholut.py nsctpl_rholut_wrap.cxx

