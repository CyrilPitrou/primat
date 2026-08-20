CC ?= cc
# -MMD/-MP emit a per-object .d file listing the headers it #includes, so a
# header edit (e.g. config.h, only #included transitively) correctly marks
# every .o that depends on it stale 
CFLAGS ?= -std=c11 -Wall -Wextra -O2 -Iinclude -MMD -MP -pthread -D_GNU_SOURCE
LDFLAGS ?= -lm -pthread
# Appended, not replacing CFLAGS: lets CI add -Werror without having to
# restate the whole flag list (and drift from it).
CFLAGS += $(CFLAGS_EXTRA)

BUILD := build
SRC := $(wildcard src/*.c)
OBJ := $(patsubst src/%.c,$(BUILD)/%.o,$(SRC))

UNIT_SRC := $(wildcard tests/unit/*.c)
UNIT_BIN := $(patsubst tests/unit/%.c,$(BUILD)/unit/%,$(UNIT_SRC))

.PHONY: all debug test debug-test bench clean

all: $(BUILD)/primat-c

debug: CFLAGS := -std=c11 -Wall -Wextra -O0 -g -fsanitize=address,undefined -Iinclude -MMD -MP -pthread -D_GNU_SOURCE
debug: LDFLAGS := -lm -pthread -fsanitize=address,undefined
debug: all

# Builds and runs the unit-test suite instrumented with ASan+UBSan (same
# flags as the `debug` target, applied to the unit binaries rather than
# $(BUILD)/primat-c). Target-specific CFLAGS/LDFLAGS propagate to the
# $(UNIT_BIN) prerequisites of `test`, so this is the cheapest memory-bug
# detector available for a C codebase with ~250 manual allocation sites and
# few NULL checks (see CI job `c-sanitizers.yml`).
debug-test: CFLAGS := -std=c11 -Wall -Wextra -O0 -g -fsanitize=address,undefined -Iinclude -MMD -MP -pthread -D_GNU_SOURCE
debug-test: LDFLAGS := -lm -pthread -fsanitize=address,undefined
debug-test: test

$(BUILD)/%.o: src/%.c
	@mkdir -p $(BUILD)
	$(CC) $(CFLAGS) -c $< -o $@

-include $(OBJ:.o=.d)

$(BUILD)/primat-c: $(OBJ)
	$(CC) $(CFLAGS) $^ -o $@ $(LDFLAGS)

$(BUILD)/unit/%: tests/unit/%.c $(filter-out $(BUILD)/main.o,$(OBJ))
	@mkdir -p $(BUILD)/unit
	$(CC) $(CFLAGS) $< $(filter-out $(BUILD)/main.o,$(OBJ)) -o $@ $(LDFLAGS)

# test_weak_rates_thermal.c #includes src/weak_rates.c directly (to reach
# its file-static L_CCRTh_compute), so it must not also link build/weak_rates.o
# (duplicate symbols) -- override the generic pattern rule above for this
# one binary.
$(BUILD)/unit/test_weak_rates_thermal: tests/unit/test_weak_rates_thermal.c $(filter-out $(BUILD)/main.o $(BUILD)/weak_rates.o,$(OBJ))
	@mkdir -p $(BUILD)/unit
	$(CC) $(CFLAGS) $< $(filter-out $(BUILD)/main.o $(BUILD)/weak_rates.o,$(OBJ)) -o $@ $(LDFLAGS)

test: $(UNIT_BIN)
	@for t in $(UNIT_BIN); do echo "== $$t =="; $$t || exit 1; done

# Times the two reference runs (small, large+amax=8) end-to-end and reports
# the ratio against the recorded Python baseline
# -- examples/baseline_timings.txt, captured once via
# `time python runfiles/primat_run.py`-equivalent calls (see that file's
# header for how/when to refresh it). `--data-dir ../primat/data` matches
# every other CPRIMAT invocation in this repo (tests, examples/*.ini);
# `mkdir -p results` only matters on a fresh checkout (cpr_nuclear_network_
# write_final_result's own mkdir_p creates it too, but cfg->output_final_
# result defaults to False via the CLI path used here -- see cli.c).
bench: $(BUILD)/primat-c
	@mkdir -p results
	@echo "== CPRIMAT benchmark (small, large+amax=8) vs. Python baseline =="
	@echo "-- small --"; \
	 start=$$(date +%s.%N); \
	 ./$(BUILD)/primat-c --data-dir ../primat/data --ini examples/run_small.ini >/dev/null; \
	 end=$$(date +%s.%N); \
	 awk -v s=$$start -v e=$$end 'BEGIN { printf "primat-c: %.3f s\n", e - s }'
	@echo "-- large, amax=8 --"; \
	 start=$$(date +%s.%N); \
	 ./$(BUILD)/primat-c --data-dir ../primat/data --ini examples/run_large_amax8.ini >/dev/null; \
	 end=$$(date +%s.%N); \
	 awk -v s=$$start -v e=$$end 'BEGIN { printf "primat-c: %.3f s\n", e - s }'
	@echo "Python baseline timings (recorded once, see examples/baseline_timings.txt):"
	@cat examples/baseline_timings.txt

# Memory-leak check: runs tests/unit/test_memory_stress.c (many repeated
# cprimat_run/cpr_mc_uncertainty/cpr_config_* init-solve-free cycles,
# including the error path) under whichever leak checker the platform
# has -- macOS's `leaks --atExit` (ships with Xcode CLT; its exit code is
# 1 iff any leak is found, 0 otherwise) or, elsewhere, valgrind
# --leak-check=full --error-exitcode=1 (same convention). A debug build
# (-O0 -g, see the `debug` target) gives valgrind/leaks better symbols and
# stack traces than the default -O2 build, so this target forces one via
# a dedicated $(BUILD)/leakcheck tree rather than reusing $(BUILD)/unit.
LEAKCHECK_BUILD := $(BUILD)/leakcheck
LEAKCHECK_CFLAGS := -std=c11 -Wall -Wextra -O0 -g -Iinclude -MMD -MP -pthread -D_GNU_SOURCE
LEAKCHECK_OBJ := $(patsubst src/%.c,$(LEAKCHECK_BUILD)/%.o,$(SRC))

$(LEAKCHECK_BUILD)/%.o: src/%.c
	@mkdir -p $(LEAKCHECK_BUILD)
	$(CC) $(LEAKCHECK_CFLAGS) -c $< -o $@

-include $(LEAKCHECK_OBJ:.o=.d)

$(LEAKCHECK_BUILD)/test_memory_stress: tests/unit/test_memory_stress.c $(filter-out $(LEAKCHECK_BUILD)/main.o,$(LEAKCHECK_OBJ))
	@mkdir -p $(LEAKCHECK_BUILD)
	$(CC) $(LEAKCHECK_CFLAGS) $< $(filter-out $(LEAKCHECK_BUILD)/main.o,$(LEAKCHECK_OBJ)) -o $@ -lm -pthread

leak-test: $(LEAKCHECK_BUILD)/test_memory_stress
	@echo "== Memory leak check (test_memory_stress) =="
	@if [ "$$(uname)" = "Darwin" ]; then \
		leaks --atExit -- $(LEAKCHECK_BUILD)/test_memory_stress; \
	elif command -v valgrind >/dev/null 2>&1; then \
		valgrind --leak-check=full --error-exitcode=1 --errors-for-leak-kinds=definite,indirect \
			$(LEAKCHECK_BUILD)/test_memory_stress; \
	else \
		echo "no leak checker available (need macOS 'leaks' or 'valgrind')"; exit 1; \
	fi

clean:
	rm -rf $(BUILD)

# ---------------------------------------------------------------------------
# Fuzzing (tests/fuzz)
#
# Apple's Command Line Tools clang instruments with -fsanitize-coverage but
# ships no libFuzzer runtime, so tests/fuzz/fuzz_engine.c supplies the loop
# (see its header). Targets keep libFuzzer's LLVMFuzzerTestOneInput signature.
# Everything under test is built with ASan+UBSan and edge instrumentation; the
# engine itself is not instrumented, so its own edges stay out of the map.
FUZZ_BUILD := $(BUILD)/fuzz
FUZZ_SAN := -fsanitize=address,undefined -fno-sanitize-recover=undefined
FUZZ_CFLAGS := -std=c11 -Wall -Wextra -O1 -g $(FUZZ_SAN) -Iinclude -Itests/fuzz \
               -MMD -MP -pthread -D_GNU_SOURCE
FUZZ_COV := -fsanitize-coverage=trace-pc-guard
FUZZ_SRC := $(wildcard tests/fuzz/fuzz_*.c)
FUZZ_TARGET_SRC := $(filter-out tests/fuzz/fuzz_engine.c,$(FUZZ_SRC))
FUZZ_BIN := $(patsubst tests/fuzz/%.c,$(FUZZ_BUILD)/%,$(FUZZ_TARGET_SRC))
FUZZ_OBJ := $(patsubst src/%.c,$(FUZZ_BUILD)/obj/%.o,$(SRC))

$(FUZZ_BUILD)/obj/%.o: src/%.c
	@mkdir -p $(FUZZ_BUILD)/obj
	$(CC) $(FUZZ_CFLAGS) $(FUZZ_COV) -c $< -o $@

$(FUZZ_BUILD)/engine.o: tests/fuzz/fuzz_engine.c
	@mkdir -p $(FUZZ_BUILD)
	$(CC) $(FUZZ_CFLAGS) -c $< -o $@

-include $(FUZZ_OBJ:.o=.d)

$(FUZZ_BUILD)/%: tests/fuzz/%.c $(FUZZ_BUILD)/engine.o $(filter-out $(FUZZ_BUILD)/obj/main.o,$(FUZZ_OBJ))
	@mkdir -p $(FUZZ_BUILD)
	$(CC) $(FUZZ_CFLAGS) $(FUZZ_COV) $< $(FUZZ_BUILD)/engine.o \
	    $(filter-out $(FUZZ_BUILD)/obj/main.o,$(FUZZ_OBJ)) -o $@ -lm -pthread $(FUZZ_SAN)

.PHONY: fuzz fuzz-build
fuzz-build: $(FUZZ_BIN)

# RUNS is per target; override on the command line (make fuzz RUNS=200000).
RUNS ?= 20000
fuzz: fuzz-build
	python3 tests/fuzz/run_fuzz.py --runs $(RUNS)

# Line/region coverage each fuzz target reached, measured by replaying its
# evolved corpus (-runs=0 runs the corpus and stops) through a second build
# instrumented with -fprofile-instr-generate instead of the sanitizers.
COV_BUILD := $(BUILD)/fuzzcov
COV_CFLAGS := -std=c11 -Wall -Wextra -O0 -g -fprofile-instr-generate \
              -fcoverage-mapping -Iinclude -Itests/fuzz -MMD -MP -pthread -D_GNU_SOURCE
COV_OBJ := $(patsubst src/%.c,$(COV_BUILD)/obj/%.o,$(SRC))
COV_BIN := $(patsubst tests/fuzz/%.c,$(COV_BUILD)/%,$(FUZZ_TARGET_SRC))

$(COV_BUILD)/obj/%.o: src/%.c
	@mkdir -p $(COV_BUILD)/obj
	$(CC) $(COV_CFLAGS) -c $< -o $@

$(COV_BUILD)/engine.o: tests/fuzz/fuzz_engine.c
	@mkdir -p $(COV_BUILD)
	$(CC) $(COV_CFLAGS) -c $< -o $@

-include $(COV_OBJ:.o=.d)

$(COV_BUILD)/%: tests/fuzz/%.c $(COV_BUILD)/engine.o $(filter-out $(COV_BUILD)/obj/main.o,$(COV_OBJ))
	@mkdir -p $(COV_BUILD)
	$(CC) $(COV_CFLAGS) $< $(COV_BUILD)/engine.o \
	    $(filter-out $(COV_BUILD)/obj/main.o,$(COV_OBJ)) -o $@ -lm -pthread \
	    -fprofile-instr-generate

.PHONY: fuzz-coverage
fuzz-coverage: $(COV_BIN)
	python3 tests/fuzz/fuzz_coverage.py
