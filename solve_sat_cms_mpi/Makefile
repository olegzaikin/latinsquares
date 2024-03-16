CPP = mpicxx
CPPFLAGS = -cxx=g++ -O3 -std=c++11

solve_sat_cms_mpi: solve_sat_cms_mpi.o
	${CPP} ${CPPFLAGS} solve_sat_cms_mpi.o -o solve_sat_cms_mpi

solve_sat_cms_mpi.o: solve_sat_cms_mpi.cpp
	${CPP} ${CPPFLAGS} solve_sat_cms_mpi.cpp -c

clean:
	rm -rf *.o
	rm solve_sat_cms_mpi
	clear
