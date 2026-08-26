/* log.c -- see log.h. */
#include "log.h"

#include <stdarg.h>
#include <stdio.h>

#if defined(_WIN32)
#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#include <io.h>       /* _isatty */
#endif

void cpr_log(const CPRConfig *cfg, const char *tag, const char *fmt, ...)
{
    if (!cfg->verbose) return;

    fprintf(stderr, "[%s-c] ", tag);
    va_list args;
    va_start(args, fmt);
    vfprintf(stderr, fmt, args);
    va_end(args);
    fputc('\n', stderr);
}

void cpr_warn(const char *fmt, ...)
{
    fputs("warning: ", stderr);
    va_list args;
    va_start(args, fmt);
    vfprintf(stderr, fmt, args);
    va_end(args);
    fputc('\n', stderr);
}

void cpr_log_raw(const CPRConfig *cfg, const char *fmt, ...)
{
    if (!cfg->verbose) return;

    va_list args;
    va_start(args, fmt);
    vfprintf(stderr, fmt, args);
    va_end(args);
    fputc('\n', stderr);
}

int cpr_console_takes_utf8(void)
{
#if defined(_WIN32)
    if (!_isatty(_fileno(stderr))) return 1;   /* redirected: raw bytes */
    return GetConsoleOutputCP() == CP_UTF8;
#else
    return 1;
#endif
}
