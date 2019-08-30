INCLUDE= -I/opt/local/include
LIBPATH= -L/opt/local/lib
CVLIBS = -lopencv_core -lopencv_highgui -lopencv_imgcodecs
FLAGS  = -O3 -Wno-unused-function
DIRBIN = ../../bin

all: convert 

convert: convert.cpp median.c
	gcc -O3 -c median.c
	g++ -c convert.cpp $(FLAGS) $(INCLUDE)
	g++ -o $(DIRBIN)/spf-convert $(FLAGS) $(LIBPATH) $(CVLIBS) median.o convert.o -lm
	rm *.o

