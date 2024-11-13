CPP = g++
CPPFLAGS = -O3 -std=c++20 -DNDEBUG

dls_main_class_enumeration: dls_main_class_enumeration.o
	${CPP} ${CPPFLAGS} dls_main_class_enumeration.o -o dls_main_class_enumeration

dls_main_class_enumeration.o: dls_main_class_enumeration.cpp
	${CPP} ${CPPFLAGS} dls_main_class_enumeration.cpp -c

clean:
	rm -rf *.o
	rm dls_main_class_enumeration
	clear
