/* fuzz.h -- the contract between a fuzz target and the engine.
 *
 * Targets define LLVMFuzzerTestOneInput (libFuzzer's signature, so they build
 * unchanged against libFuzzer where its runtime is available) and get a
 * private scratch directory to materialise the input as a file, which is what
 * every parser under test actually consumes.
 */
#ifndef CPRIMAT_FUZZ_H
#define CPRIMAT_FUZZ_H

#include <stddef.h>
#include <stdint.h>

int LLVMFuzzerTestOneInput(const uint8_t *data, size_t size);

/* Called once before the first input. Every target defines it -- Mach-O has
 * no portable weak-undefined symbol, so a no-op body is cheaper than the
 * linker flags that would make it optional. */
int LLVMFuzzerInitialize(int *argc, char ***argv);

/* Per-process scratch directory, created by the engine and removed at exit. */
const char *fuzz_tmpdir(void);

/* Write `size` bytes to <tmpdir>/<relpath>, creating parent directories.
 * Returns a pointer to a static buffer holding the full path. */
const char *fuzz_write_file(const char *relpath, const uint8_t *data, size_t size);

/* Data root the targets read shipped tables from (CPRIMAT_FUZZ_DATA_DIR,
 * default ../primat/data). */
const char *fuzz_data_dir(void);

#endif /* CPRIMAT_FUZZ_H */
