/* test_xalloc.c -- exercises the success paths of the checked-allocation
 * helpers cpr_xmalloc/cpr_xcalloc/cpr_xrealloc (see ../../include/xalloc.h).
 *
 * The failure path (log + exit on genuine OOM) cannot be exercised in-process
 * without forking or interposing malloc, and its only observable effect is a
 * clean exit(1) with a diagnostic -- so this test pins the parts that run on
 * every solve: that a normal request returns usable, correctly-sized memory,
 * that cpr_xcalloc zero-initialises, that cpr_xrealloc grows a block while
 * preserving its prefix, that cpr_xrealloc(NULL, n) behaves like a fresh
 * allocation, and that a zero-size request does not spuriously abort. Run
 * under ASan/UBSan in CI, this also confirms the returned blocks are valid,
 * writable, and freeable across their full extent. */
#include "xalloc.h"

#include <stdio.h>
#include <stdlib.h>

static int failures = 0;

#define CHECK(cond, msg) do { \
        if (!(cond)) { printf("FAIL: %s\n", msg); failures++; } \
        else printf("ok: %s\n", msg); \
    } while (0)

int main(void)
{
    /* cpr_xmalloc: non-NULL, and the whole block is writable (ASan checks bounds). */
    size_t n = 256;
    double *v = CPR_XMALLOC(n * sizeof(double));
    CHECK(v != NULL, "cpr_xmalloc returns non-NULL");
    for (size_t i = 0; i < n; i++) v[i] = (double)i;
    CHECK(v[0] == 0.0 && v[n - 1] == (double)(n - 1), "cpr_xmalloc block is fully writable");
    free(v);

    /* cpr_xcalloc: non-NULL and zero-initialised across the full extent. */
    size_t *c = CPR_XCALLOC(n, sizeof(size_t));
    CHECK(c != NULL, "cpr_xcalloc returns non-NULL");
    int all_zero = 1;
    for (size_t i = 0; i < n; i++) if (c[i] != 0) all_zero = 0;
    CHECK(all_zero, "cpr_xcalloc zero-initialises the block");

    /* cpr_xrealloc: grows the block and preserves the existing prefix. */
    for (size_t i = 0; i < n; i++) c[i] = i + 1;
    c = CPR_XREALLOC(c, 2 * n * sizeof(size_t));
    CHECK(c != NULL, "cpr_xrealloc returns non-NULL");
    int prefix_ok = 1;
    for (size_t i = 0; i < n; i++) if (c[i] != i + 1) prefix_ok = 0;
    CHECK(prefix_ok, "cpr_xrealloc preserves the original prefix");
    for (size_t i = n; i < 2 * n; i++) c[i] = 0; /* touch the grown tail */
    free(c);

    /* cpr_xrealloc(NULL, n) is a fresh allocation (like malloc). */
    char *fresh = CPR_XREALLOC(NULL, 64);
    CHECK(fresh != NULL, "cpr_xrealloc(NULL, n) allocates like malloc");
    fresh[0] = 'x'; fresh[63] = 'y';
    free(fresh);

    /* Zero-size request must not abort; its (possibly-NULL) result is freeable. */
    void *z = CPR_XMALLOC(0);
    free(z);
    CHECK(1, "cpr_xmalloc(0) does not abort");

    if (failures) { printf("\n%d FAILURE(S)\n", failures); return 1; }
    printf("\nAll xalloc tests passed.\n");
    return 0;
}
