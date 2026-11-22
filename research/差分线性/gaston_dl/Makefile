CXX ?= g++
CC ?= cc
CPPFLAGS += -Iinclude
CXXFLAGS ?= -O3 -std=c++17 -Wall -Wextra -Wpedantic
LDLIBS += -pthread
REF_DIR ?= ../../gaston_ref/gaston/implementation/reference

.PHONY: all test test-six-six test-pq-kernel-six test-datf test-datf-fixed-tail test-datf-pools test-linear-two-pool test-datf-quadratic test-datf-paper test-fixed-tail-extend test-hnw-extension sanitize
all: build/gaston_dl build/fixed_tail_verify build/conditional_trail_verify

build:
	mkdir -p build

build/gaston_dl: src/main.cpp include/gaston.hpp include/estimator.hpp include/experiment.hpp | build
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@ $(LDLIBS)

build/fixed_tail_verify: src/fixed_tail_verify.cpp include/gaston.hpp include/estimator.hpp include/experiment.hpp | build
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@ $(LDLIBS)

build/conditional_trail_verify: src/conditional_trail_verify.cpp include/gaston.hpp include/estimator.hpp include/experiment.hpp | build
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@ $(LDLIBS)

build/datf_middle_cli: src/datf_middle_cli.cpp include/datf_middle.hpp include/gaston.hpp | build
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@

build/datf_pool_match: src/datf_pool_match.cpp include/datf_middle.hpp include/gaston.hpp | build
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@

build/datf_paper_shard: src/datf_paper_shard.cpp include/datf_middle.hpp include/gaston.hpp | build
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@

build/datf_linear_two_shard: src/datf_linear_two_shard.cpp include/datf_middle.hpp include/gaston.hpp | build
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@

build/hnw_endpoint: src/hnw_endpoint.cpp include/gaston.hpp include/estimator.hpp | build
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@ -lgmpxx -lgmp

test-hnw-extension: build/hnw_endpoint
	./build/hnw_endpoint --selftest
	python3 tests/test_hnw_extension.py

build/datf_pool_verify: src/datf_pool_verify.cpp include/experiment.hpp include/gaston.hpp | build
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@ $(LDLIBS)

build/xof_interface_baseline: src/xof_interface_baseline.cpp include/gaston.hpp build/reference.o
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -I$(REF_DIR) $< build/reference.o -o $@ $(LDLIBS)

.PHONY: test-xof-interface
test-xof-interface: build/xof_interface_baseline
	python3 tests/test_xof_interface.py

build/fixed_tail_extend: src/fixed_tail_extend.cpp include/gaston.hpp | build
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@

test-fixed-tail-extend: build/fixed_tail_extend
	./build/fixed_tail_extend --self-test
	python3 tests/test_fixed_tail_extend.py

build/tests_datf_middle: tests/test_datf_middle.cpp include/datf_middle.hpp include/gaston.hpp | build
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@

build/reference.o: $(REF_DIR)/gaston.c $(REF_DIR)/gaston.h | build
	$(CC) -O2 -I$(REF_DIR) -c $< -o $@

build/tests: tests/test.cpp include/gaston.hpp include/estimator.hpp include/experiment.hpp build/reference.o
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -I$(REF_DIR) $< build/reference.o -o $@ $(LDLIBS)

test: build/tests
	./build/tests
	python3 tests/test_fixed_tail.py

test-six-six: build/conditional_trail_verify
	python3 tests/test_six_six.py

test-pq-kernel-six:
	python3 tests/test_pq_kernel_six.py

test-datf:
	python3 tests/test_datf.py

test-datf-fixed-tail:
	python3 tests/test_datf_fixed_tail.py

test-datf-pools: build/datf_middle_cli build/tests_datf_middle
	./build/tests_datf_middle
	python3 tests/test_datf_pools.py

test-linear-two-pool: build/datf_pool_match build/datf_middle_cli
	python3 tests/test_linear_two_pool.py

test-datf-quadratic:
	python3 tests/test_datf_quadratic.py

test-datf-paper:
	python3 tests/test_datf_paper_full.py

sanitize: build/reference.o
	$(CXX) $(CPPFLAGS) -I$(REF_DIR) -O1 -g -std=c++17 -fsanitize=address,undefined -fno-omit-frame-pointer tests/test.cpp build/reference.o -o build/tests-sanitize $(LDLIBS)
	./build/tests-sanitize
