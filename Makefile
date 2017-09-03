CXX = g++
CXXFLAGS= -Wall -g #-std=c++11 #-Wall #-g
LINKPATH= -I./samtools-0.1.19 -L./samtools-0.1.19
LINKFLAGS = -lbam -lz -lm -lpthread 
DEBUG=
OBJECTS = stats.o subexon-graph.o 

all: subexon-info combine-subexons classes grader

subexon-info: subexon-info.o $(OBJECTS)
	$(CXX) -o $@ $(LINKPATH) $(CXXFLAGS) $(OBJECTS) subexon-info.o $(LINKFLAGS)

combine-subexons: combine-subexons.o $(OBJECTS)
	$(CXX) -o $@ $(LINKPATH) $(CXXFLAGS) $(OBJECTS) combine-subexons.o $(LINKFLAGS)

classes: classes.o constraints.o transcript-decider.o $(OBJECTS)
	$(CXX) -o $@ $(LINKPATH) $(CXXFLAGS) $(OBJECTS) constraints.o transcript-decider.o classes.o $(LINKFLAGS)

grader: grader.o
	$(CXX) -o $@ $(LINKPATH) $(CXXFLAGS) $(OBJECTS) grader.o $(LINKFLAGS)
	

subexon-info.o: SubexonInfo.cpp alignments.hpp blocks.hpp support.hpp defs.h stats.hpp
	$(CXX) -c -o $@ $(LINKPATH) $(CXXFLAGS) $< $(LINKFLAGS)
combine-subexons.o: CombineSubexons.cpp alignments.hpp blocks.hpp support.hpp defs.h stats.hpp SubexonGraph.hpp
	$(CXX) -c -o $@ $(LINKPATH) $(CXXFLAGS) $< $(LINKFLAGS)
stats.o: stats.cpp stats.hpp
	$(CXX) -c -o $@ $(LINKPATH) $(CXXFLAGS) $< $(LINKFLAGS)
subexon-graph.o: SubexonGraph.cpp SubexonGraph.hpp
	$(CXX) -c -o $@ $(LINKPATH) $(CXXFLAGS) $< $(LINKFLAGS)
constraints.o: Constraints.cpp Constraints.hpp alignments.hpp
	$(CXX) -c -o $@ $(LINKPATH) $(CXXFLAGS) $< $(LINKFLAGS)
transcript-decider.o: TranscriptDecider.cpp TranscriptDecider.hpp Constraints.hpp BitTable.hpp
	$(CXX) -c -o $@ $(LINKPATH) $(CXXFLAGS) $< $(LINKFLAGS)
classes.o: classes.cpp SubexonGraph.hpp SubexonCorrelation.hpp BitTable.hpp Constraints.hpp alignments.hpp TranscriptDecider.hpp
	$(CXX) -c -o $@ $(LINKPATH) $(CXXFLAGS) $< $(LINKFLAGS)
grader.o: grader.cpp
	$(CXX) -c -o $@ $(LINKPATH) $(CXXFLAGS) $< $(LINKFLAGS)

clean:
	rm -f *.o *.gch subexon-info combine-subexons 
