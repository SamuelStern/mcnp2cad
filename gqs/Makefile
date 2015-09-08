

include ${CGM_MAKE}

#include ${MOAB_MAKE}

all: make_gq

test: test_hyperbolic_curves

gq.o: 
	g++ -g -c gq.cpp -o gq.o $(CGM_INCLUDES) $(CGM_LIBS_LINK)

make_gq: gq.o ProgOptions.o
	g++ -g make_gq.cpp gq.o ProgOptions.o -o make_gq  $(CGM_INCLUDES) -I../ $(CGM_LIBS_LINK)

ProgOptions.o:
	g++ -c ../ProgOptions.cpp ../ProgOptions.hpp

test_hyperbolic_curves: gq.o ProgOptions.o
	g++ -g test_hyperbolic_curves.cpp gq.o ProgOptions.o -o test_hyperbolic_curves $(CGM_INCLUDES) $(CGM_LIBS_LINK)

clean:
	rm gq.o make_gq test_hyperbolic_curves
