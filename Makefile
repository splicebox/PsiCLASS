CXX = g++
CXXFLAGS= #-Wall #-g
LINKPATH= -I./samtools-0.1.19 -L./samtools-0.1.19
LINKFLAGS = -lbam -lz -lm -lpthread 
DEBUG=
OBJECTS = gamma.o

all: subexon-info

subexon-info: main.o $(OBJECTS)
	$(CXX) -o $@ $(LINKPATH) $(CXXFLAGS) $(OBJECTS) main.o $(LINKFLAGS)

main.o: main.cpp alignments.hpp blocks.hpp support.hpp defs.h gamma.o
gamma.o: gamma.cpp gamma.hpp

clean:
	rm -f *.o *.gch subexon-info 
