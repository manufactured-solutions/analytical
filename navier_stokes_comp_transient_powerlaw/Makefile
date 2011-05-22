all: writeup.pdf check

writeup.pdf: writeup.tex soln.py forcing.py nsctpl_fwd.hpp nsctpl.hpp
	md5sum soln.py forcing.py nsctpl_fwd.hpp nsctpl.hpp > CHECKSUMS
	pdflatex writeup
	pdflatex writeup

test.hpp: test.py soln.py forcing.py
	test.py | tee test.hpp

check: check.cpp nsctpl_fwd.hpp nsctpl.hpp test.hpp
	g++ -g -O0 -fno-unsafe-math-optimizations -Wall -pedantic check.cpp -o check

module: nsctpl_fwd.hpp nsctpl.hpp nsctpl.i
	swig -v -Wall -c++ -python nsctpl.i
	python setup.py build_ext --inplace

clean:
	rm -fv *.o check test.hpp soln.c soln.h *.aux *.log *.out writeup.pdf

