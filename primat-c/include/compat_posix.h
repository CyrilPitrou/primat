#ifndef CPRIMAT_COMPAT_POSIX_H
#define CPRIMAT_COMPAT_POSIX_H

/* =========================================================================
 * compat_posix.h -- POSIX/Windows portability shim (non-threading part).
 *
 * The C backend uses a handful of POSIX facilities (the <unistd.h>/
 * <sys/stat.h> headers and a few libc functions) that either do not exist
 * or are spelled differently under MSVC.  Rather than sprinkle #ifdef _WIN32
 * blocks through every source file, each affected .c file includes this
 * single header in place of <unistd.h>/<sys/stat.h>:
 *
 *   - POSIX (Linux/macOS): we simply pull in the real <unistd.h> and
 *     <sys/stat.h>, so behaviour is byte-for-byte unchanged.
 *   - Windows/MSVC: we include the CRT headers that declare the equivalents
 *     and map the POSIX spellings onto them via macros, so the rest of the
 *     sources stay platform-agnostic.
 *
 * The threading primitives (pthreads) are shimmed separately in
 * compat_thread.h, which mc.c includes.
 * ========================================================================= */

#if defined(_WIN32)

#include <direct.h>   /* _mkdir, _getcwd                                    */
#include <io.h>       /* _access (not currently used, kept for parity)      */
#include <process.h>  /* _getpid                                            */
#include <string.h>   /* strtok_s (MSVC's re-entrant strtok, == strtok_r)   */
#include <sys/stat.h> /* struct stat / stat(): MSVC ships this header too   */

/* strtok_s has the exact signature of POSIX strtok_r:
 *   char *strtok_s(char *str, const char *delim, char **context);          */
#ifndef strtok_r
#define strtok_r strtok_s
#endif

/* Case-insensitive string compares: POSIX str[n]casecmp -> MSVC _str[n]icmp
 * (identical signatures). */
#ifndef strcasecmp
#define strcasecmp _stricmp
#endif
#ifndef strncasecmp
#define strncasecmp _strnicmp
#endif

/* POSIX mkdir(path, mode) -> MSVC _mkdir(path): the CRT variant takes no
 * permission-mode argument (Windows ACLs are inherited), so we drop it.    */
#ifndef mkdir
#define mkdir(path, mode) _mkdir(path)
#endif

/* getcwd / getpid live under leading-underscore names in the MSVC CRT.     */
#ifndef getcwd
#define getcwd _getcwd
#endif
#ifndef getpid
#define getpid _getpid
#endif

#else /* POSIX (Linux, macOS, ...) */

#include <unistd.h>
#include <sys/stat.h>

#endif /* _WIN32 */

#endif /* CPRIMAT_COMPAT_POSIX_H */
