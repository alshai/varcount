CXX=g++
CXX_FLAGS=-std=c++11 -O3 -Wall -Wextra -Ihtslib

LIBS=-lhts
CPP=$(wildcard *.cpp)
OBJ=$(CPP:%.cpp=%.o)

all: varcount compare_vcfs

varcount: varcount.o hts_util.hpp
	$(CXX) $(CXX_FLAGS) $^ -o $@ $(LIBS)

compare_vcfs: compare_vcfs.o hts_util.hpp
	$(CXX) $(CXX_FLAGS) $^ -o $@ $(LIBS)

test: test.o hts_util.hpp
	$(CXX) $(CXX_FLAGS) $^ -o $@ $(LIBS)

%.o: %.cpp
	$(CXX) $(CXX_FLAGS) -c $< -o $@ $(LIBS)

.PHONY: clean

clean: 
	-rm varcount test $(OBJ)
