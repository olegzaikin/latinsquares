CPP = g++
CPPFLAGS = -O3 -std=c++20

find_pseudotriple_dls_dmd_style: find_pseudotriple_dls_dmd_style.o dlx_orth.o
	${CPP} ${CPPFLAGS} find_pseudotriple_dls_dmd_style.o dlx_orth.o -o find_pseudotriple_dls_dmd_style

find_pseudotriple_dls_dmd_style.o: find_pseudotriple_dls_dmd_style.cpp
	${CPP} ${CPPFLAGS} find_pseudotriple_dls_dmd_style.cpp -c

dlx_orth.o: ../dlx_mols/dlx_orth.cpp
	${CPP} ${CPPFLAGS}  ../dlx_mols/dlx_orth.cpp -c

clean:
	rm -rf *.o
	rm find_pseudotriple_dls_dmd_style
	clear
