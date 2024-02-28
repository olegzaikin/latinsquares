CPP = g++
CPPFLAGS = -O3 -std=c++20 -fopenmp

solve_sat_cms: solve_sat_cms.o
	${CPP} ${CPPFLAGS} solve_sat_cms.o -o solve_sat_cms

solve_sat_cms.o: solve_sat_cms.cpp
	${CPP} ${CPPFLAGS} solve_sat_cms.cpp -c

clean:
	rm -rf *.o
	rm solve_sat_cms
	clear
