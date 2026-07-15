/* xalloc.c -- see xalloc.h for the full rationale.
 *
 * These are deliberately tiny: allocate via libc, and on the (rare, genuinely
 * unrecoverable) NULL return, emit a precise diagnostic and stop. The
 * diagnostic goes to stderr and stdout is flushed first, so that any partial
 * verbose (`cpr_log`) output already produced is not lost or interleaved after
 * the error line. */
#include "xalloc.h"

#include <stdio.h>
#include <stdlib.h>

/* Shared failure path: report the byte count and source location, then exit.
 * Marked so it never returns; kept out of the header because it is an internal
 * implementation detail of the three wrappers. */
static void cpr_alloc_fail(size_t size, const char *file, int line)
{
    fflush(stdout);
    fprintf(stderr, "primat: out of memory (%zu bytes) at %s:%d\n",
            size, file, line);
    exit(EXIT_FAILURE);
}

void *cpr_xmalloc(size_t size, const char *file, int line)
{
    void *p = malloc(size);
    /* A zero-size request may legitimately return NULL on some libcs; only
     * treat NULL as failure when we actually asked for memory. */
    if (p == NULL && size != 0)
        cpr_alloc_fail(size, file, line);
    return p;
}

void *cpr_xcalloc(size_t nmemb, size_t size, const char *file, int line)
{
    void *p = calloc(nmemb, size);
    if (p == NULL && nmemb != 0 && size != 0)
        cpr_alloc_fail(nmemb * size, file, line);
    return p;
}

void *cpr_xrealloc(void *ptr, size_t size, const char *file, int line)
{
    void *p = realloc(ptr, size);
    if (p == NULL && size != 0)
        cpr_alloc_fail(size, file, line);
    return p;
}
