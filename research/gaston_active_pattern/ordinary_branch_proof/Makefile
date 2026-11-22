CXX ?= g++
SRC ?= enumerate_effects.cpp
TARGET ?= $(basename $(SRC))
CXXFLAGS ?= -O3 -march=native -std=c++17 -Wall -Wextra
GUROBI_HOME ?= /opt/gurobi1302/linux64
.PHONY: all
all: $(TARGET)
ifneq ($(filter ordinary_milp effect_cover,$(TARGET)),)
CPPFLAGS += -I$(GUROBI_HOME)/include
LDFLAGS += -L$(GUROBI_HOME)/lib -Wl,-rpath,$(GUROBI_HOME)/lib
LDLIBS += -lgurobi_c++ -lgurobi130 -pthread
endif
$(TARGET): $(SRC) ../direct_enumeration/gaston_common.hpp activity_matching.hpp
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $< $(LDFLAGS) $(LDLIBS) -o $@
%: %.cpp
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $< $(LDFLAGS) $(LDLIBS) -o $@
