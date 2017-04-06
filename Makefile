CXX = g++
CXXFLAGS= -std=c++0x -Wall #-std=c++11 #-Wall #-g
LINKPATH= -I./samtools-0.1.19 -L./samtools-0.1.19
LINKFLAGS = -lbam -lz -lm -lpthread 
DEBUG=
OBJECTS = gamma.o

all: subexon-info combineSubexons

subexon-info: main.o $(OBJECTS)
	$(CXX) -o $@ $(LINKPATH) $(CXXFLAGS) $(OBJECTS) main.o $(LINKFLAGS)

combineSubexons: CombineSubexons.cpp alignments.hpp $(OBJECTS)
	$(CXX) -o $@ $(LINKPATH) $(CXXFLAGS) $(OBJECTS) CombineSubexons.cpp $(LINKFLAGS)

main.o: main.cpp alignments.hpp blocks.hpp support.hpp defs.h gamma.cpp gamma.hpp
gamma.o: gamma.cpp gamma.hpp

clean:
	rm -f *.o *.gch subexon-info combineSubexons 
