all: main pde

main: main.C dualnumber.h numberarray.h compare_types.h
	g++ -o main main.C 

pde: pde.C dualnumberarray.h dualnumber.h numberarray.h compare_types.h
	g++ -o pde pde.C -I$(MASA_INC) -L$(MASA_LIB) -lmasa

clean: 
	rm pde main
