all: main pde

main: main.C dualnumber.h numberarray.h
	g++ -o main main.C 

pde: pde.C dualnumberarray.h
	g++ -o pde pde.C -I$(MASA_INC) -L$(MASA_LIB) -lmasa

dualnumber.h: compare_types.h
dualnumber.h: raw_type.h

number_array.h: compare_types.h
number_array.h: raw_type.h

dualnumberarray.h: dualnumber.h
dualnumberarray.h: number_array.h

clean: 
	rm pde main