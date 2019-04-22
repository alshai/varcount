CXX=g++
CXX_FLAGS=--std=c++11 -O3 -Wall -Wextra

LIBS=-lhts
CPP=$(wildcard *.cpp)
OBJ=$(CPP:%.cpp=%.o)
DEP=$(OBJ:%.o=%.d)

all: varcount

varcount: varcount.o hts_util.o
	$(CXX) $(CXX_FLAGS) $^ -o $@ $(LIBS)

-include $(DEP)

%.o: %.cpp
	$(CXX) $(CXX_FLAGS) -MMD -c $< -o $@ $(LIBS)

.PHONY: clean

clean: 
	-rm varcount $(OBJ) $(DEP)
