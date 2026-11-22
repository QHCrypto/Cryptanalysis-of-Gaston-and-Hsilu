CXX ?= g++
GUROBI_HOME ?= /opt/gurobi1302/linux64
SRC ?= hsilu_3round_milp.cpp
TARGET ?= $(basename $(notdir $(SRC)))

CXXFLAGS ?= -O2 -std=c++17 -Wall -Wextra -Wpedantic
CPPFLAGS += -I$(GUROBI_HOME)/include
LDFLAGS += -L$(GUROBI_HOME)/lib -Wl,-rpath,$(GUROBI_HOME)/lib
LDLIBS += -lgurobi_c++ -lgurobi130 -pthread
REFERENCE = ../hsilu_ref/hsilu.cpp

.PHONY: all clean
all: $(TARGET)

$(TARGET): $(SRC) $(REFERENCE) ../hsilu_ref/hsilu.h
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(SRC) $(REFERENCE) $(LDFLAGS) $(LDLIBS) -o $@

clean:
	$(RM) $(TARGET)
