CPP = g++
CPPFLAGS = -O3 -std=c++20

sat_enumeration_brown_dls: sat_enumeration_brown_dls.o dlx_orth.o
	${CPP} ${CPPFLAGS} sat_enumeration_brown_dls.o dlx_orth.o -o sat_enumeration_brown_dls

sat_enumeration_brown_dls.o: sat_enumeration_brown_dls.cpp
	${CPP} ${CPPFLAGS} sat_enumeration_brown_dls.cpp -c

dlx_orth.o: ../dlx_mols/dlx_orth.cpp
	${CPP} ${CPPFLAGS}  ../dlx_mols/dlx_orth.cpp -c

clean:
	rm -rf *.o
	rm sat_enumeration_brown_dls
	clear
