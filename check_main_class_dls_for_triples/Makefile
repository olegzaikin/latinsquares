CPP = g++
CPPFLAGS = -O3 -std=c++20

check_main_class_dls_for_triples: check_main_class_dls_for_triples.o
	${CPP} ${CPPFLAGS} check_main_class_dls_for_triples.o -o check_main_class_dls_for_triples

check_main_class_dls_for_triples.o: check_main_class_dls_for_triples.cpp
	${CPP} ${CPPFLAGS} check_main_class_dls_for_triples.cpp -c

clean:
	rm -rf *.o
	rm check_main_class_dls_for_triples
	clear
