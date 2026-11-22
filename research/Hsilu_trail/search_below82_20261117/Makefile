CXX ?= g++
CXXFLAGS ?= -O3 -march=native -std=c++17 -Wall -Wextra

search: search.cpp common.hpp
	$(CXX) $(CXXFLAGS) search.cpp -o search

clean:
	rm -f search
