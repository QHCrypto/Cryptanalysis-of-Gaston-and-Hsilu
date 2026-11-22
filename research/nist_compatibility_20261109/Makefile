CXX ?= g++
CXXFLAGS ?= -O3 -std=c++17 -Wall -Wextra -Wpedantic

.PHONY: all verify
all: native.so

native.so: native.cpp ../hsilu_ref/hsilu.cpp ../hsilu_ref/hsilu.h ../gaston分析/three_round_algebra.hpp
	$(CXX) $(CXXFLAGS) -shared -fPIC native.cpp ../hsilu_ref/hsilu.cpp -o $@

verify: native.so
	python3 audit.py
