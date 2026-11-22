CXX ?= g++
GUROBI_HOME ?= /opt/gurobi1302/linux64
SRC ?= hsilu_2round_milp.cpp
TARGET ?= $(basename $(notdir $(SRC)))

CXXFLAGS ?= -O3 -std=c++17 -Wall -Wextra -Wpedantic
CPPFLAGS += -I$(GUROBI_HOME)/include
LDFLAGS += -L$(GUROBI_HOME)/lib -Wl,-rpath,$(GUROBI_HOME)/lib
LDLIBS += -lgurobi_c++ -lgurobi130 -pthread

.PHONY: all clean check
all: $(TARGET)

$(TARGET): $(SRC)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< $(LDFLAGS) $(LDLIBS) -o $@

%: %.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< $(LDFLAGS) $(LDLIBS) -o $@

check: hsilu_2round_milp
	./hsilu_2round_milp --check-reference hsilu_2round_best_attack.txt --digest-bits 64 --output reproduced_reference --threads 4
	python3 validate_milp_path.py reproduced_reference.txt

clean:
	$(RM) $(TARGET)
