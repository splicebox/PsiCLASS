CXX = g++
CXXFLAGS= -Wall #-std=c++11 #-Wall #-g
LINKPATH= -I./samtools-0.1.19 -L./samtools-0.1.19
LINKFLAGS = -lbam -lz -lm -lpthread 
DEBUG=
OBJECTS = stats.o

all: subexon-info combine-subexons

subexon-info: subexon-info.o stats.o
	$(CXX) -o $@ $(LINKPATH) $(CXXFLAGS) $(OBJECTS) subexon-info.o $(LINKFLAGS)

combine-subexons: combine-subexons.o $(OBJECTS)
	$(CXX) -o $@ $(LINKPATH) $(CXXFLAGS) $(OBJECTS) combine-subexons.o $(LINKFLAGS)

subexon-info.o: SubexonInfo.cpp alignments.hpp blocks.hpp support.hpp defs.h stats.cpp stats.hpp
	$(CXX) -c -o $@ $(LINKPATH) $(CXXFLAGS) $< $(LINKFLAGS)
combine-subexons.o: CombineSubexons.cpp alignments.hpp blocks.hpp support.hpp defs.h stats.cpp stats.hpp
	$(CXX) -c -o $@ $(LINKPATH) $(CXXFLAGS) $< $(LINKFLAGS)
stats.o: stats.cpp stats.hpp
	$(CXX) -c -o $@ $(LINKPATH) $(CXXFLAGS) $< $(LINKFLAGS)

clean:
	rm -f *.o *.gch subexon-info combineSubexons 
