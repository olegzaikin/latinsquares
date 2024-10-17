CPP = g++
CPPFLAGS = -O3 -std=c++20

dlx_mols: dlx_mols.o dlx_orth.o
	${CPP} ${CPPFLAGS} dlx_mols.o dlx_orth.o -o dlx_mols

dlx_mols.o: dlx_mols.cpp
	${CPP} ${CPPFLAGS} dlx_mols.cpp -c

dlx_orth.o: dlx_orth.cpp
	${CPP} ${CPPFLAGS} dlx_orth.cpp -c

clean:
	rm -rf *.o
	rm dlx_mols
	clear
