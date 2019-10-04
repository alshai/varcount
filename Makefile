CXX=g++
WARNINGS=-Wall -Wextra -Wno-char-subscripts \
         -Wpointer-arith -Wwrite-strings -Wdisabled-optimization \
         -Wformat -Wcast-align -Wno-unused-function -Wno-unused-parameter \
         -pedantic -Wunused-variable -Wno-attributes -Wno-pedantic  -Wno-ignored-attributes
CXX_FLAGS=-std=c++11 -O3 -Wall -Wextra -Ihtslib $(WARNINGS)

LIBS=-lz -lbz2 -llzma -lcurl
CPP=$(wildcard *.cpp)
OBJ=$(CPP:%.cpp=%.o)

all: varcount next_gt test unit

htslib/libhts.a:
	cd htslib && make libhts.a

varcount: varcount.o hts_util.hpp
	$(CXX) $(CXX_FLAGS) $^ htslib/libhts.a -o $@ $(LIBS)

compare_vcfs: compare_vcfs.o hts_util.hpp
	$(CXX) $(CXX_FLAGS) $^ htslib/libhts.a -o $@ $(LIBS)

vcf_score: vcf_score.o hts_util.hpp
	$(CXX) $(CXX_FLAGS) $^ htslib/libhts.a -o $@ $(LIBS)

next_gt: next_gt.o hts_util.hpp
	$(CXX) $(CXX_FLAGS) $^ htslib/libhts.a -o $@ $(LIBS)

test: test.o hts_util.hpp
	$(CXX) $(CXX_FLAGS) $^ htslib/libhts.a -o $@ $(LIBS)

unit: testutil.o
	$(CXX) $< -o $@ $(CXX_FLAGS) $(LIBS)

testutil: testutil.o
	$(CXX) $< -o $@ $(CXX_FLAGS) $(LIBS)

%.o: %.cpp
	$(CXX) $(CXX_FLAGS) -c $< -o $@ $(LIBS)

.PHONY: clean

clean:
	-rm varcount test $(OBJ)
