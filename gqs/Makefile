

include ${CGM_MAKE}

include ${MOAB_MAKE}

all: gq

gq.o: 
	g++ -g -c gq.cpp -o gq.o $(CGM_INCLUDES) $(CGM_LIBS_LINK)
make_gq: gq.o
	g++ -g make_gq.cpp gq.hpp gq.o -o make_gq  $(CGM_INCLUDES) $(MOAB_INCLUDES) \
	$(CGM_LIBS_LINK) $(MOAB_LIBS_LINK)
test: test_hyperbolic_curves

test_hyperbolic_curves:
	g++ -g gq.hpp gq.cpp test_hyperbolic_curves.cpp -o test_hyperbolic_curves $(CGM_INCLUDES) $(MOAB_INCLUDES) \
	$(CGM_LIBS_LINK) $(MOAB_LIBS_LINK)
clean:
	rm gq.o make_gq test_hyperbolic_curves
