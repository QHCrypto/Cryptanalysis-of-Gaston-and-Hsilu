CXX ?= g++
CC ?= gcc
GASTON_REF := ../../gaston_ref/gaston/implementation/reference
CXXFLAGS ?= -O3 -std=c++17 -Wall -Wextra -Wpedantic -flto
CFLAGS ?= -O3 -std=c11 -Wall -Wextra -flto

all: certificate kernel_probe orbit_coefficients exterior_zdd

gaston_reference.o: $(GASTON_REF)/gaston.c $(GASTON_REF)/gaston.h
	$(CC) $(CFLAGS) -I$(GASTON_REF) -c $< -o $@

certificate: certificate.cpp gaston_reference.o
	$(CXX) $(CXXFLAGS) -I$(GASTON_REF) $^ -o $@

kernel_probe: kernel_probe.cpp gaston_reference.o
	$(CXX) $(CXXFLAGS) -I$(GASTON_REF) $^ -o $@

orbit_coefficients: orbit_coefficients.cpp
	$(CXX) $(CXXFLAGS) -I$(GASTON_REF) $< -o $@

exterior_zdd: exterior_zdd.cpp
	$(CXX) $(CXXFLAGS) $< -o $@

.PHONY: all clean
clean:
	rm -f certificate kernel_probe orbit_coefficients exterior_zdd gaston_reference.o
