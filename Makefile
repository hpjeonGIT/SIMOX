#
# Makefile for SIMOX serial version
# 
# Byoungseon Jeon
# School of Engineering and Applied Sciences
# Harvard University
# July 2010
#
.SUFFIXES: .o .f90
#
F90 =  ifort #gfortran
#FLAGS =  -O3 -m64 -ipo -no-prec-div #-heap-arrays  
#FLAGS =  -march=native -ffast-math -funroll-loops -O3 
#FLAGS =  -g -Wall -O3
#FLAGS =  -xHost -O3 -ipo -no-prec-div  -pc80 -pad -ip
FLAGS = -xHOST  -O3 -ipo -no-prec-div  # -heap-arrays
F77FLAG = ${FLAGS} #-O3 -m64 -ipo -no-prec-div #-heap-arrays  # 
#F77FLAG = -march=native -ffast-math -funroll-loops -O3 -fomit-frame-pointer -I..
OBJ = datafmt.o main.o force.o vverlet.o ctip.o util.o eam.o cellsort.o \
      fft235.o kernel.o mfft235.o zfft3d.o
#LIB = -L/opt/local/atlas/lib/ -llapack -lf77blas -lcblas -latlas
TARGET = simox_cellsort
${TARGET}:${OBJ}
	${F90} ${FLAGS} -o ${TARGET} ${OBJ} ${LIB} 

.f90.o:
	${F90} ${FLAGS} -c $< ${INC}
.f.o:
	${F90} ${F77FLAG} -c $<
clean:
	rm -rf *.mod *.o *.f90~ core ${TARGET}
