/* xalloc.h -- checked heap allocation for the primat-c backend.
 *
 * Motivation
 * ----------
 * The solver allocates a few hundred small/medium buffers per run (ODE work
 * vectors, spline coefficient tables, the loaded reaction network, the
 * time-evolution result arrays, ...). The overwhelming majority of these sit
 * in pure numerical kernels (ode_bdf/ode_rk/spline/linalg/network_builder/
 * plasma/qed_pressure/weak_rates result assembly/...) that have *no* graceful
 * recovery path: if the allocation fails, the very next line dereferences the
 * pointer and the process crashes with an uninformative SIGSEGV -- or, worse,
 * silently corrupts memory. For a scientific CLI, aborting on genuine
 * out-of-memory is entirely acceptable; silently dereferencing NULL is not.
 *
 * These helpers centralise the "allocate, and if it fails log where and die"
 * contract in one place, so those unrecoverable call sites become a single
 * self-documenting token (CPR_XMALLOC(...) instead of malloc(...)) and every
 * OOM produces a precise `file:line` diagnostic instead of a bare crash.
 *
 * When NOT to use these
 * ---------------------
 * Call sites that are *meant* to degrade gracefully -- above all the cache
 * writers, which warn and continue when a cache file cannot be produced --
 * must keep their explicit `if (p == NULL) { warn; ... }` handling and use the
 * plain libc malloc/calloc/realloc. `cpr_x*` is for the "there is nothing
 * sensible to do but stop" majority only.
 *
 * Note on the Python extension
 * ----------------------------
 * primat-c is also linked into the CPython extension (primat._primat_c). A
 * failed allocation there will terminate the host interpreter via exit(3)
 * rather than raising a Python exception. That is a deliberate trade-off: a
 * true OOM deep inside the BDF integrator has no clean unwinding path across
 * the C ABI anyway, and a labelled exit is strictly better than the NULL
 * dereference these call sites exhibited before. Functions that already
 * propagate a recoverable error to Python via a `char **errmsg` out-parameter
 * keep doing so (they are in the "graceful" bucket above and are left alone).
 *
 * Usage
 * -----
 *     double *v = CPR_XMALLOC(n * sizeof(double));   // never NULL on return
 *     size_t *c = CPR_XCALLOC(n, sizeof(size_t));    // zero-initialised
 *     buf       = CPR_XREALLOC(buf, new_cap);        // grow an existing block
 *     ... use v/c/buf ...
 *     free(v); free(c); free(buf);
 *
 * The macros capture __FILE__/__LINE__ automatically; the underlying
 * functions can also be called directly if a call site needs to pass an
 * explicit location.
 */
#ifndef CPRIMAT_XALLOC_H
#define CPRIMAT_XALLOC_H

#include <stddef.h>

/* Allocate `size` bytes. On success returns a non-NULL pointer to
 * uninitialised memory (like malloc). On failure prints
 * "primat: out of memory (<size> bytes) at <file>:<line>" to stderr and
 * exit(EXIT_FAILURE)s -- it never returns NULL. A zero `size` is forwarded to
 * malloc unchanged (its implementation-defined result is treated as success,
 * matching libc behaviour; callers that care must guard `size > 0` themselves). */
void *cpr_xmalloc(size_t size, const char *file, int line);

/* Allocate and zero `nmemb * size` bytes (like calloc), with the same
 * never-returns-NULL / log-and-exit-on-failure contract as cpr_xmalloc. */
void *cpr_xcalloc(size_t nmemb, size_t size, const char *file, int line);

/* Resize the block `ptr` to `size` bytes (like realloc), with the same
 * log-and-exit-on-failure contract. `ptr == NULL` behaves like cpr_xmalloc.
 * On failure the original block is left untouched (it is leaked, which is
 * irrelevant since the process exits immediately). */
void *cpr_xrealloc(void *ptr, size_t size, const char *file, int line);

/* Convenience macros that stamp in the current source location. Prefer these
 * over calling the functions directly so the OOM diagnostic points at the
 * real call site. */
#define CPR_XMALLOC(size)        cpr_xmalloc((size), __FILE__, __LINE__)
#define CPR_XCALLOC(nmemb, size) cpr_xcalloc((nmemb), (size), __FILE__, __LINE__)
#define CPR_XREALLOC(ptr, size)  cpr_xrealloc((ptr), (size), __FILE__, __LINE__)

#endif /* CPRIMAT_XALLOC_H */
