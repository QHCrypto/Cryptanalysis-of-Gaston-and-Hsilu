GUROBI_HOME ?= /opt/gurobi1302/linux64
CXX ?= g++
CXXFLAGS ?= -std=c++17 -O2 -Wall -Wextra -Wpedantic
CPPFLAGS += -I$(GUROBI_HOME)/include
LDFLAGS += -L$(GUROBI_HOME)/lib -Wl,-rpath,$(GUROBI_HOME)/lib
LDLIBS += -lgurobi_c++ -lgurobi130 -pthread

SRC ?= original_milp.cpp
TARGET ?= $(basename $(notdir $(SRC)))

.PHONY: all both all_models check_three_round check_value clean
all: $(TARGET)
both: original_milp improved_milp
all_models: original_milp improved_milp value_milp three_round_milp
check_three_round: three_round_algebra_test three_round_analysis_test
	./three_round_algebra_test
	./three_round_analysis_test

original_milp improved_milp: model_common.hpp
three_round_milp: three_round_analysis.hpp three_round_algebra.hpp
value_milp: value_algebra.hpp three_round_algebra.hpp
check_value: value_algebra_test
	./value_algebra_test

value_algebra_test: value_algebra_test.cpp value_algebra.hpp three_round_algebra.hpp
	$(CXX) $(CXXFLAGS) $< -o $@

%: %.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $@.cpp $(LDFLAGS) $(LDLIBS) -o $@

three_round_algebra_test: three_round_algebra_test.cpp three_round_algebra.hpp
	$(CXX) $(CXXFLAGS) $< -o $@

three_round_analysis_test: three_round_analysis_test.cpp three_round_analysis.hpp three_round_algebra.hpp
	$(CXX) $(CXXFLAGS) $< -o $@

clean:
	$(RM) original_milp improved_milp value_milp value_algebra_test three_round_milp three_round_algebra_test three_round_analysis_test
