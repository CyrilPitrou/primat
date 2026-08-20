/* A target that crashes on a known input, so the harness can be shown to see
 * a crash at all. "no crashes in N executions" means nothing from a pipeline
 * that cannot report one; run_fuzz.py --selftest drives this target and fails
 * unless an artifact reproducing the fault comes back. */
#include "fuzz.h"

#include <stdlib.h>
#include <string.h>

int LLVMFuzzerInitialize(int *argc, char ***argv)
{
    (void)argc; (void)argv;
    return 0;
}

int LLVMFuzzerTestOneInput(const uint8_t *data, size_t size)
{
    /* A heap overflow ASan reports, on a 4-byte token a mutator reaches in a
     * few thousand execs from the seeded corpus. The read goes through a
     * volatile sink: an allocate-copy-free triple whose result is unused is
     * dead code, and clang deletes it outright at -O1 -- which is what the
     * first version of this target did, so it never crashed. */
    if (size >= 4 && memcmp(data, "BOOM", 4) == 0) {
        char *p = malloc(size);
        memcpy(p, data, size);
        volatile char sink = p[size];   /* one past the end */
        (void)sink;
        free(p);
    }
    return 0;
}
