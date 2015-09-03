

include ${CGM_MAKE}

include ${MOAB_MAKE}

all: gq hyperbola

gq:
	g++ -g gq.hpp  gq.cpp -o gq $(CGM_INCLUDES) $(MOAB_INCLUDES) $(CGM_LIBS_LINK) $(MOAB_LIBS_LINK)

hyperbola:
	g++ -g gq.hpp hyperbola.cpp -o hyperbola $(CGM_INCLUDES) $(MOAB_INCLUDES) $(CGM_LIBS_LINK) $(MOAB_LIBS_LINK)

clean:
	rm gq hyperbola
