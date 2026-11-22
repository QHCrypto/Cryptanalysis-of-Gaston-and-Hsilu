CXX ?= g++
CXXFLAGS ?= -O3 -std=c++17 -Wall -Wextra -Wpedantic
REF_DIR ?= ../../hsilu_ref

.PHONY: all test
all: build/hsilu_dl

build:
	mkdir -p build

build/hsilu_dl: src/main.cpp include/hsilu.hpp | build
	$(CXX) $(CXXFLAGS) -Iinclude -I$(REF_DIR) $< $(REF_DIR)/hsilu.cpp -pthread -o $@

test: build/hsilu_dl
	./build/hsilu_dl selftest
	python3 tests/test_core.py
	python3 tests/test_prefix.py

build/hnw_endpoint: src/hnw_endpoint.cpp include/hnw_adapter.hpp include/hsilu.hpp | build
	$(CXX) $(CXXFLAGS) -Iinclude $< -lgmpxx -lgmp -o $@

build/reference_verify: src/reference_verify.cpp | build
	$(CXX) $(CXXFLAGS) -I$(REF_DIR) $< $(REF_DIR)/hsilu.cpp -o $@
