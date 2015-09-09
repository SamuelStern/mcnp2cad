

include ${CGM_MAKE}

#include ${MOAB_MAKE}
SHELL = /bin/bash -o pipefail

all: make_gq

tests: test_hyperbolic_curves test_hyperbolic_curves_in_plane test_types

check: tests
	./run_tests | tee test-log.txt
gq.o: 
	g++ -g -c gq.cpp -o gq.o $(CGM_INCLUDES) $(CGM_LIBS_LINK)

make_gq: gq.o ProgOptions.o
	g++ -g make_gq.cpp gq.o ProgOptions.o -o make_gq  $(CGM_INCLUDES) -I../ $(CGM_LIBS_LINK)

ProgOptions.o:
	g++ -c ../ProgOptions.cpp ../ProgOptions.hpp

test_hyperbolic_curves: gq.o ProgOptions.o
	g++ -g test_hyperbolic_curves.cpp gq.o ProgOptions.o -o test_hyperbolic_curves $(CGM_INCLUDES) $(CGM_LIBS_LINK)

test_hyperbolic_curves_in_plane: gq.o ProgOptions.o
	g++ -g test_hyperbolic_curves_in_plane.cpp gq.o ProgOptions.o -o test_hyperbolic_curves_in_plane $(CGM_INCLUDES) $(CGM_LIBS_LINK)

test_types: gq.o ProgOptions.o
	g++ -g test_types.cpp gq.o ProgOptions.o -o test_types $(CGM_INCLUDES) $(CGM_LIBS_LINK)

clean:
	rm gq.o make_gq test_hyperbolic_curves test_hyperbolic_curves_in_plane test_types  *.sat *.txt
