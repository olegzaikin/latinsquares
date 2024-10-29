CPP = g++
CPPFLAGS = -O3 -std=c++20

enumerate_brown_dls: enumerate_brown_dls.o dlx_orth.o
	${CPP} ${CPPFLAGS} enumerate_brown_dls.o dlx_orth.o -o enumerate_brown_dls

enumerate_brown_dls.o: enumerate_brown_dls.cpp
	${CPP} ${CPPFLAGS} enumerate_brown_dls.cpp -c

dlx_orth.o: ../dlx_mols/dlx_orth.cpp
	${CPP} ${CPPFLAGS}  ../dlx_mols/dlx_orth.cpp -c

clean:
	rm -rf *.o
	rm enumerate_brown_dls
	clear
