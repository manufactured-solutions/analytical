CXX ?= g++

all: main pde

main: main.C dualnumber.h numberarray.h shadownumber.h dualshadowarray.h compare_types.h raw_type.h
	$(CXX) -Wall -o main main.C 

pde: pde.C dualnumberarray.h dualnumber.h numberarray.h compare_types.h raw_type.h shadownumber.h dualshadowarray.h
	$(CXX) -Wall -o pde pde.C -I$(MASA_INC) -L$(MASA_LIB) -lmasa

clean: 
	rm pde main
