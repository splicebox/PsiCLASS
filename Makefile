CXX = g++
CXXFLAGS= -Wall #-std=c++11 #-Wall #-g
LINKPATH= -I./samtools-0.1.19 -L./samtools-0.1.19
LINKFLAGS = -lbam -lz -lm -lpthread 
DEBUG=
OBJECTS = stats.o

all: subexon-info combineSubexons

subexon-info: main.o $(OBJECTS)
	$(CXX) -o $@ $(LINKPATH) $(CXXFLAGS) $(OBJECTS) main.o $(LINKFLAGS)

combineSubexons: CombineSubexons.cpp blocks.hpp alignments.hpp $(OBJECTS)
	$(CXX) -o $@ $(LINKPATH) $(CXXFLAGS) $(OBJECTS) CombineSubexons.cpp $(LINKFLAGS)

main.o: main.cpp alignments.hpp blocks.hpp support.hpp defs.h stats.cpp stats.hpp
stats.o: stats.cpp stats.hpp

clean:
	rm -f *.o *.gch subexon-info combineSubexons 
