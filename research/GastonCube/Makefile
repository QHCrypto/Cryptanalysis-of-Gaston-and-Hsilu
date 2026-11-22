CXX ?= g++
GUROBI_HOME ?= /opt/gurobi1302/linux64
SRC ?= gaston_cube.cpp
TARGET ?= $(basename $(SRC))
CXXFLAGS ?= -O3 -std=c++17 -Wall -Wextra -Wpedantic
CPPFLAGS += -I$(GUROBI_HOME)/include
LDFLAGS += -L$(GUROBI_HOME)/lib -Wl,-rpath,$(GUROBI_HOME)/lib
LDLIBS += -lgurobi_c++ -lgurobi130 -pthread

.PHONY: all test clean
all: $(TARGET)
$(TARGET): $(SRC) gaston_mp.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< $(LDFLAGS) $(LDLIBS) -o $@
%: %.cpp gaston_mp.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< $(LDFLAGS) $(LDLIBS) -o $@
test: all gaston_parity_cube
	python3 tests/validate.py
	python3 tests/validate_parity.py
clean:
	rm -f $(TARGET)
