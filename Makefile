CXX=g++
CXX_FLAGS=--std=c++11 -O3 -Wall -Wextra

LIBS=-lhts
CPP=$(wildcard *.cpp)
OBJ=$(CPP:%.cpp=%.o)

all: varcount 

varcount: varcount.o hts_util.hpp
	$(CXX) $(CXX_FLAGS) $^ -o $@ $(LIBS)

test: test.o hts_util.hpp
	$(CXX) $(CXX_FLAGS) $^ -o $@ $(LIBS)

%.o: %.cpp
	$(CXX) $(CXX_FLAGS) -c $< -o $@ $(LIBS)

.PHONY: clean

clean: 
	-rm varcount test $(OBJ)
