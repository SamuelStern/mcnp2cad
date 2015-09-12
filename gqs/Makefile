
SHELL = /bin/bash -o pipefail
include ${CGM_MAKE}
include ${MOAB_MAKE}

all: make_gq

tests: test_hyperbolic_curves test_hyperbolic_curves_in_plane test_types

check: tests
	./run_tests | tee test-log.txt

make_gq: ProgOptions.o
	g++ -g make_gq.cpp gq.cpp mcnp2cad_funcs.cpp ProgOptions.o -o make_gq  $(CGM_INCLUDES) $(MOAB_INCLUDES) -I../ $(CGM_LIBS_LINK) $(MOAB_LIBS_LINK)

ProgOptions.o:
	g++ -c ../ProgOptions.cpp ../ProgOptions.hpp

test_hyperbolic_curves: ProgOptions.o
	g++ -g test_hyperbolic_curves.cpp gq.cpp mcnp2cad_funcs.cpp  ProgOptions.o -o test_hyperbolic_curves $(CGM_INCLUDES) $(MOAB_INCLUDES) $(CGM_LIBS_LINK) $(MOAB_LIBS_LINK)

test_hyperbolic_curves_in_plane:  ProgOptions.o
	g++ -g test_hyperbolic_curves_in_plane.cpp gq.cpp mcnp2cad_funcs.cpp  ProgOptions.o -o test_hyperbolic_curves_in_plane $(CGM_INCLUDES) $(MOAB_INCLUDES) $(CGM_LIBS_LINK) $(MOAB_LIBS_LINK)

test_types:  ProgOptions.o
	g++ -g test_types.cpp gq.cpp mcnp2cad_funcs.cpp ProgOptions.o -o test_types $(CGM_INCLUDES) $(MOAB_INCLUDES) $(CGM_LIBS_LINK) $(MOAB_LIBS_LINK)

clean:
	rm make_gq test_hyperbolic_curves test_hyperbolic_curves_in_plane test_types  *.sat *.txt
