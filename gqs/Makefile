

include ${CGM_MAKE}

include ${MOAB_MAKE}

all: gq

gq:
	g++ -g gq.cpp -o gq $(CGM_INCLUDES) $(MOAB_INCLUDES) $(CGM_LIBS_LINK) $(MOAB_LIBS_LINK)
clean:
	rm gq
