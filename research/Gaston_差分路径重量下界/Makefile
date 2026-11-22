CXX ?= g++
GUROBI_HOME ?= /opt/gurobi1302/linux64
SRC ?= gaston_3round_lossy.cpp
TARGET ?= $(basename $(SRC))
CXXFLAGS ?= -O3 -std=c++17 -Wall -Wextra -Wpedantic
CPPFLAGS += -I$(GUROBI_HOME)/include
LDFLAGS += -L$(GUROBI_HOME)/lib -Wl,-rpath,$(GUROBI_HOME)/lib
LDLIBS += -lgurobi_c++ -lgurobi130 -pthread

.PHONY: all clean support-proof-binaries
all: $(TARGET)
support-proof-binaries: src/support_rank_probe src/verify_support_ranks
src/support_rank_probe: src/support_rank_probe.cpp
	$(CXX) $(CXXFLAGS) $< -o $@
src/verify_support_ranks: src/verify_support_ranks.cpp
	$(CXX) $(CXXFLAGS) $< -o $@
$(TARGET): $(SRC) src/activity_profile_data.hpp src/kernel_type_data.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< $(LDFLAGS) $(LDLIBS) -o $@
%: %.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< $(LDFLAGS) $(LDLIBS) -o $@
clean:
	rm -f $(TARGET)
