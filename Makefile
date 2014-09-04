
#MAIN = test
#MAIN = main_u1
#MAIN = main_dirac
#MAIN = main_u1real
#MAIN = main_u1compl
#MAIN = main_u1exact
#MAIN = main_u1real_scan12
#MAIN = main_u1compl_scan12
#MAIN = main_u1kagome_scan12
#MAIN = main_u1compl_scan13
#MAIN = main_twoflavor
#MAIN = main_u1hyb2
#MAIN = main_he_col
#MAIN = main_he
MAIN = main_he_kag2

CLASSES = helperfunctions.cpp isingstate.cpp subspace.cpp lattice.cpp wavefunction.cpp vmc.cpp matrix.cpp mysql_wrapper.cpp spinone.cpp u1.cpp twoflavor.cpp u1two.cpp u1hybrid.cpp u1hybtwo.cpp paired2k.cpp u1kagome.cpp u1real.cpp huseelser.cpp he_two.cpp
#u1dirac.cpp

#LFLAGS = lib/libtmglib.a lib/liblapack.a lib/libblas.a lib/libmysqlcppconn.so -lm -lgfortran
LFLAGS = lib/libpfapack.a -llapack -lblas -lmysqlcppconn -lm -lgfortran

#use e-fence or valgrind to debug memory problems
#option: -lefence or lib/libefence.a

LIB = /usr/lib64

GCC = g++

all:
	@echo "Main is: " ${MAIN}
	rm -f main*.o
#	${GCC} -c -Wall -O3 -std=gnu++11 ${CLASSES} ${MAIN}.cpp -I.
	${GCC} -c -Wall -O3 -Wno-vla -std=gnu++11 ${CLASSES} ${MAIN}.cpp -I.
	${GCC} -Wall -O3 -L${LIB} *.o -I. -o bin/${MAIN} ${LFLAGS}

fast:
	@echo "Main is: " ${MAIN}
	rm -f main*.o
	${GCC} -c -Wall -O0 -Wno-vla -std=gnu++11 ${CLASSES} ${MAIN}.cpp -I.
	${GCC} -Wall -O0 -L${LIB} *.o -I. -o bin/${MAIN} ${LFLAGS}

link:
	${GCC} -g0 -O3 -L${LIB} *.o -I. -o bin/${MAIN} ${LFLAGS}

debug:
	${GCC} -g -ggdb -Wall ${CLASSES} ${MAIN}.cpp -I. -c
	${GCC} -g -ggdb -Wall *.o -I. -o bin/${MAIN}_debug ${LFLAGS}
	gdb ${MAIN}_debug

clean:
	rm -f *.o

#performance profiling with gnu-profiler
prof:
	rm -f main*.o
	${GCC} -Wall -pg -g -O3 ${CLASSES} ${MAIN}.cpp -I. -c
	${GCC} -Wall -pg -g -O3 -L${LIB} *.o -I. -o bin/${MAIN}_debug ${LFLAGS}

valgrind:
	rm -f main*.o
	${GCC} -c -Wall -g -O1 ${CLASSES} ${MAIN}.cpp -I.
	${GCC} -Wall -g -O1 -L${LIB} *.o -I. -o bin/${MAIN}_debug ${LFLAGS}

