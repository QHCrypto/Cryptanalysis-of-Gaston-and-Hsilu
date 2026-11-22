CXX ?= g++
SRC ?= gaston_enumerate.cpp
TARGET ?= $(basename $(SRC))
CXXFLAGS ?= -O3 -march=native -std=c++17 -Wall -Wextra -Wpedantic

.PHONY: all check check-parallel clean
all: $(TARGET)

ifeq ($(filter gaston_parallel gaston_out6_nonzero,$(TARGET)),)
$(TARGET): $(SRC) gaston_common.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< $(LDFLAGS) $(LDLIBS) -o $@
endif

gaston_parallel: gaston_parallel.cpp gaston_common.hpp
	$(CXX) $(CPPFLAGS) -Ithird_party/OpenCL $(CXXFLAGS) -fopenmp -pthread $< $(LDFLAGS) -l:libOpenCL.so.1 $(LDLIBS) -o $@

gaston_out6_nonzero: gaston_out6_nonzero.cpp gaston_common.hpp
	$(CXX) $(CPPFLAGS) -Ithird_party/OpenCL $(CXXFLAGS) $< $(LDFLAGS) -l:libOpenCL.so.1 $(LDLIBS) -o $@

%: %.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< $(LDFLAGS) $(LDLIBS) -o $@

check: $(TARGET)
	python3 verify.py ./$(TARGET)

check-parallel: gaston_parallel
	python3 verify_parallel.py

clean:
	rm -f $(TARGET)
