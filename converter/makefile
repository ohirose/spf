INCLUDE= -I/opt/local/include
LIBPATH= -L/opt/local/lib
CVLIBS = -lopencv_core -lopencv_highgui
FLAGS  = -Wall -O3 -Wno-unused-function 

all: convert 

convert: convert.c median.c
	gcc -o convert $(FLAGS) $(INCLUDE) $(LIBPATH) $(CVLIBS) convert.c median.c


