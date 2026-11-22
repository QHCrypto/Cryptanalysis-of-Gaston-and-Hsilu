CXX ?= g++
GUROBI_HOME ?= /opt/gurobi1302/linux64

SRC ?= gaston_active_pattern.cpp
TARGET ?= $(basename $(SRC))

CXXFLAGS ?= -O3 -std=c++17 -Wall -Wextra -Wpedantic
CPPFLAGS += -I$(GUROBI_HOME)/include
LDFLAGS += -L$(GUROBI_HOME)/lib -Wl,-rpath,$(GUROBI_HOME)/lib
LDLIBS += -lgurobi_c++ -lgurobi130 -pthread

.PHONY: all check clean

all: $(TARGET)

$(TARGET): $(SRC)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< $(LDFLAGS) $(LDLIBS) -o $@

%: %.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< $(LDFLAGS) $(LDLIBS) -o $@

check: $(TARGET)
	python3 verify.py ./$(TARGET)

clean:
	rm -f $(TARGET)
