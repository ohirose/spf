C = gcc -Wall
INCLUDE = -I/opt/local/include/
LIBPATH = -L/opt/local/lib/
FLAG    = -Wno-deprecated-declarations
COMMON  = ../common/
DIRBIN  = ../../bin

LIBS = -lm -lglut -lGL -lGLU

Viewer : view.c 
	$(C) view.c $(COMMON)io.c $(COMMON)util.c $(FLAG) $(INCLUDE) $(LIBPATH) $(LIBS) -O3 -o $(DIRBIN)/spf-view

clean : 
	rm -f *~ src/*~ obj/*.o view





