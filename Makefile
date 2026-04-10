ifndef CC
override CC = gcc
endif

ifndef CXX
override CXX = g++
endif

ifndef ENABLE_TEST_ALL
	DEFINED_FLAGS =
else
	DEFINED_FLAGS = -DENABLE_TEST_ALL
endif

ifndef CROSS_COMPILE
    processor := $(shell uname -m)
	ARCH_CFLAGS = -maes -mpclmul -mssse3 -msse4.2
else # CROSS_COMPILE was set
    CC = $(CROSS_COMPILE)gcc
    CXX = $(CROSS_COMPILE)g++
    CXXFLAGS += -static
    LDFLAGS += -static

    check_riscv := $(shell echo | $(CROSS_COMPILE)cpp -dM - | grep " __riscv_xlen " | cut -c22-)
    uname_result := $(shell uname -m)
	ifeq ($(check_riscv),64)
		processor = rv64
    else ifeq ($(uname_result),rv64imafdc)
		processor = rv64
    else ifeq ($(check_riscv),32)
		processor = rv32
    else ifeq ($(uname_result),rv32i)
		processor = rv32
	else
		$(error Unsupported cross-compiler)
	endif

	ifeq ($(processor),$(filter $(processor),i386 x86_64))
		ARCH_CFLAGS = -maes -mpclmul -mssse3 -msse4.2
	else
		ARCH_CFLAGS = -march=$(processor)gcv_zba
	endif

	ifeq ($(SIMULATOR_TYPE), qemu)
		SIMULATOR += qemu-riscv64
		SIMULATOR_FLAGS = -cpu $(processor),v=true,zba=true,vlen=128
	else
		SIMULATOR = spike
		SIMULATOR_FLAGS = --isa=$(processor)gcv_zba
		PROXY_KERNEL = pk
	endif
endif

CXXFLAGS += -Wall -Wcast-qual -I. $(ARCH_CFLAGS)
LDFLAGS	+= -lm
OBJS = \
    tests/binding.o \
    tests/common.o \
    tests/debug_tools.o \
    tests/impl.o \
    tests/main.o
deps := $(OBJS:%.o=%.o.d)

.SUFFIXES: .o .cpp
.cpp.o:
	$(CXX) -o $@ $(CXXFLAGS) $(DEFINED_FLAGS) -c -MMD -MF $@.d $<

EXEC = tests/main
COMPAT_TEST = tests/compat_test

$(EXEC): $(OBJS)
	$(CXX) $(LDFLAGS) -o $@ $^

$(COMPAT_TEST): tests/compat_test.cpp sse2rvv.h
	$(CXX) -o $@ $(CXXFLAGS) $(LDFLAGS) tests/compat_test.cpp

test: tests/main
ifeq ($(processor),$(filter $(processor),rv32 rv64))
	$(CC) $(ARCH_CFLAGS) -c sse2rvv.h
endif
	$(SIMULATOR) $(SIMULATOR_FLAGS) $(PROXY_KERNEL) $^

compat-test: $(COMPAT_TEST)
ifeq ($(processor),$(filter $(processor),rv32 rv64))
	$(SIMULATOR) $(SIMULATOR_FLAGS) $(PROXY_KERNEL) $^
else
	$^
endif

build-test: tests/main
ifeq ($(processor),$(filter $(processor),rv32 rv64))
	$(CC) $(ARCH_CFLAGS) -c sse2rvv.h
endif

format:
	@echo "Formatting files with clang-format.."
	@if ! hash clang-format; then echo "clang-format is required to indent"; fi
	clang-format -i sse2rvv.h tests/*.cpp tests/*.h

.PHONY: clean check format compat-test

clean:
	$(RM) $(OBJS) $(EXEC) $(COMPAT_TEST) $(deps) sse2rvv.h.gch

clean-all: clean
	$(RM) *.log

-include $(deps)
