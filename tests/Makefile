

all: check


check:
	${CXX} -g ${CXXFLAGS} ${IGEOM_CPPFLAGS} -I../ ../geometry.o ../MCNPInput.o ../volumes.o test_gqs.cpp -o test_gqs ${LDFLAGS}
