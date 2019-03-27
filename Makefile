CXX=g++
CXX_FLAGS=--std=c++11 -O3

LIBS=-lhts

all: varcount

varcount: varcount.cpp mdparse.hpp varcount.hpp
	$(CXX) $(CXX_FLAGS) -o $@ $< $(LIBS)
