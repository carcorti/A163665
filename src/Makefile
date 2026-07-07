CC ?= gcc
CSTD := -std=c17
WARN := -Wall -Wextra -Wshadow -Wformat=2 -Wconversion -Wdouble-promotion -Wpedantic -Wstrict-prototypes -Wmissing-prototypes
OPT ?= -O2
OPENMP ?= -fopenmp
override CFLAGS := $(CSTD) $(OPT) $(WARN) $(OPENMP) $(EXTRA_CFLAGS)
override LDFLAGS := $(OPENMP) $(EXTRA_LDFLAGS)
override LDLIBS := -lm $(EXTRA_LDLIBS)

TARGET := a163665
SRC := a163665.c
SOURCE_PATH := $(abspath $(SRC))
SOURCE_SHA256 := $(shell sha256sum $(SRC) | awk '{print $$1}')
override CPPFLAGS := -DA163665_SOURCE_PATH=\"$(SOURCE_PATH)\" -DA163665_BUILD_SOURCE_SHA256=\"$(SOURCE_SHA256)\" $(EXTRA_CPPFLAGS)
SMOKE_PREFIX_DEFAULT := ../runs/a163665_smoke_$(shell date +%s%N)
SMOKE_PREFIX ?= $(SMOKE_PREFIX_DEFAULT)
THREADS ?= 4

.PHONY: all clean test selftest smoke sanitize

all: $(TARGET)

$(TARGET): $(SRC)
	$(CC) $(CPPFLAGS) $(CFLAGS) $< -o $@ $(LDFLAGS) $(LDLIBS)

selftest: $(TARGET)
	./$(TARGET) selftest

smoke: $(TARGET)
	mkdir -p ../runs
	./$(TARGET) segment --lo 1 --hi 100000 --threads 1 --out $(SMOKE_PREFIX)_t1
	./$(TARGET) segment --lo 1 --hi 100000 --threads $(THREADS) --out $(SMOKE_PREFIX)_tN
	cmp $(SMOKE_PREFIX)_t1.hits.tsv $(SMOKE_PREFIX)_tN.hits.tsv
	./$(TARGET) segment --lo 1 --hi 50000 --threads 1 --out $(SMOKE_PREFIX)_sum_a
	./$(TARGET) segment --lo 50001 --hi 100000 --threads $(THREADS) --out $(SMOKE_PREFIX)_sum_b
	./$(TARGET) summarize $(SMOKE_PREFIX)_sum_a.manifest.tsv $(SMOKE_PREFIX)_sum_b.manifest.tsv

test: selftest smoke

sanitize:
	$(CC) $(CPPFLAGS) $(CSTD) -O1 -g3 $(WARN) $(OPENMP) -fsanitize=address,undefined $(SRC) -o $(TARGET)_san $(LDFLAGS) -fsanitize=address,undefined $(LDLIBS)
	env -u LD_PRELOAD ASAN_OPTIONS=detect_leaks=0 ./$(TARGET)_san selftest
	env -u LD_PRELOAD ASAN_OPTIONS=detect_leaks=0 ./$(TARGET)_san segment --lo 1 --hi 20000 --threads 2 --out /tmp/a163665_san_$$(date +%s%N)

clean:
	rm -f $(TARGET) $(TARGET)_san
