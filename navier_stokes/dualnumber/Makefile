CXX ?= icpc

all: main pde pde_debug

main: main.C dualnumber.h numberarray.h shadownumber.h dualshadowarray.h compare_types.h raw_type.h
	$(CXX) $(CXXFLAGS) -Wall -o main main.C 

pde: pde.C dualnumberarray.h dualnumber.h numberarray.h compare_types.h raw_type.h shadownumber.h dualshadowarray.h
	$(CXX) $(CXXFLAGS) -O0 -Wall -o pde pde.C -I$(MASA_INC) -L$(MASA_LIB) -lmasa

pde_debug: pde_debug.C dualnumberarray.h dualnumber.h numberarray.h compare_types.h raw_type.h shadownumber.h dualshadowarray.h
	$(CXX) $(CXXFLAGS) -O0 -Wall -o pde_debug pde_debug.C -I$(MASA_INC) -L$(MASA_LIB) -lmasa

clean: 
	rm pde pde_debug main
