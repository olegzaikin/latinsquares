CPP = g++
CPPFLAGS = -O3 -DNDEBUG -std=c++20 -fopenmp

sat_enumeration_brown_dls: sat_enumeration_brown_dls.o
	${CPP} ${CPPFLAGS} sat_enumeration_brown_dls.o -o sat_enumeration_brown_dls

sat_enumeration_brown_dls.o: sat_enumeration_brown_dls.cpp
	${CPP} ${CPPFLAGS} sat_enumeration_brown_dls.cpp -c

clean:
	rm -rf *.o
	rm sat_enumeration_brown_dls
	clear
