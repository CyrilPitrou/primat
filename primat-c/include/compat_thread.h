#ifndef CPRIMAT_COMPAT_THREAD_H
#define CPRIMAT_COMPAT_THREAD_H

/* =========================================================================
 * compat_thread.h -- POSIX/Windows threading portability shim.
 *
 * The Monte-Carlo driver (mc.c) parallelises sample evaluation across worker
 * threads using the POSIX pthreads subset {mutex, create, join} plus
 * sysconf(_SC_NPROCESSORS_ONLN) for CPU autodetection.  MSVC has none of
 * these, so on Windows we provide a minimal drop-in implementation backed by
 * Win32 primitives; on POSIX we forward to the real <pthread.h>/<unistd.h>.
 *
 * Only the facilities the callers actually use are shimmed -- this is not a
 * general pthreads emulation.  Two translation units include this header
 * (mc.c and primat/_primat_c_src/_wrapper.c), each using a different subset,
 * so the definitions are `static inline` -- unused ones in either unit must
 * not trigger -Wunused-function.
 * ========================================================================= */

#if defined(_WIN32)

#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#endif
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>
#include <process.h>  /* _beginthreadex */
#include <stdint.h>   /* uintptr_t      */
#include <stdlib.h>   /* malloc/free    */

/* --- Opaque handle types -------------------------------------------------
 * A thread is a Win32 HANDLE; a mutex is a CRITICAL_SECTION (lighter than a
 * kernel mutex and sufficient for intra-process locking, which is all mc.c
 * needs).                                                                   */
typedef HANDLE pthread_t;
typedef CRITICAL_SECTION pthread_mutex_t;

/* --- Mutex -------------------------------------------------------------- */
static inline int pthread_mutex_init(pthread_mutex_t *m, void *attr) {
    (void)attr; /* mc.c always passes NULL (default attributes) */
    InitializeCriticalSection(m);
    return 0;
}
static inline int pthread_mutex_destroy(pthread_mutex_t *m) {
    DeleteCriticalSection(m);
    return 0;
}
static inline int pthread_mutex_lock(pthread_mutex_t *m) {
    EnterCriticalSection(m);
    return 0;
}
static inline int pthread_mutex_unlock(pthread_mutex_t *m) {
    LeaveCriticalSection(m);
    return 0;
}

/* --- Thread create/join --------------------------------------------------
 * A pthread start routine has signature `void *(*)(void *)`, whereas
 * _beginthreadex expects `unsigned __stdcall (*)(void *)`.  We bridge the two
 * with a small heap-allocated trampoline record carrying the real function
 * pointer and argument.  _beginthreadex (not raw CreateThread) is used so the
 * per-thread C runtime state is initialised/torn down correctly -- the worker
 * bodies call malloc/free and other CRT functions.                          */
typedef void *(*cpr_pthread_fn)(void *);
struct cpr_pthread_start {
    cpr_pthread_fn fn;
    void *arg;
};

static inline unsigned __stdcall cpr_pthread_trampoline(void *p) {
    struct cpr_pthread_start s = *(struct cpr_pthread_start *)p;
    free(p);
    s.fn(s.arg); /* return value is discarded; mc.c passes NULL to join */
    return 0;
}

static inline int pthread_create(pthread_t *t, void *attr,
                          cpr_pthread_fn fn, void *arg) {
    (void)attr; /* mc.c always passes NULL (default attributes) */
    struct cpr_pthread_start *s =
        (struct cpr_pthread_start *)malloc(sizeof *s);
    if (!s) {
        return -1;
    }
    s->fn = fn;
    s->arg = arg;
    uintptr_t h = _beginthreadex(NULL, 0, cpr_pthread_trampoline, s, 0, NULL);
    if (h == 0) {
        free(s);
        return -1;
    }
    *t = (HANDLE)h;
    return 0;
}

static inline int pthread_join(pthread_t t, void **retval) {
    (void)retval; /* mc.c never inspects worker return values */
    WaitForSingleObject(t, INFINITE);
    CloseHandle(t);
    return 0;
}

/* --- CPU count -----------------------------------------------------------
 * sysconf(_SC_NPROCESSORS_ONLN) -> number of logical processors.            */
#ifndef _SC_NPROCESSORS_ONLN
#define _SC_NPROCESSORS_ONLN 1
#endif
static inline long sysconf(int name) {
    (void)name; /* mc.c only ever asks for _SC_NPROCESSORS_ONLN */
    SYSTEM_INFO si;
    GetSystemInfo(&si);
    return (long)si.dwNumberOfProcessors;
}

/* --- usleep --------------------------------------------------------------
 * _wrapper.c's Ctrl-C poll loop sleeps in microseconds; Win32 Sleep() takes
 * milliseconds.  Round up so a sub-millisecond request still yields the CPU. */
static inline int usleep(unsigned long usec) {
    Sleep((DWORD)((usec + 999UL) / 1000UL));
    return 0;
}

#else /* POSIX (Linux, macOS, ...) */

#include <pthread.h>
#include <unistd.h>

#endif /* _WIN32 */

#endif /* CPRIMAT_COMPAT_THREAD_H */
