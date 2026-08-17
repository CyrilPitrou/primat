/* log.h -- verbose-mode console logging, the C-side counterpart of
 * primat's scattered `if cfg.verbose: print(...)` lines (see documentation's
 * "Unify and complete --verbose output" effort). Centralised here so every
 * call site is a single line and the `[<tag>-c] ` prefix convention (mirrors
 * Python's `[<tag>-py] ` tags) stays consistent.
 */
#ifndef CPRIMAT_LOG_H
#define CPRIMAT_LOG_H

#include "config.h"

/* No-op unless cfg->verbose. Prints "[<tag>-c] " followed by the formatted
 * message and a trailing newline to stdout. `tag` should be one of
 * "init"/"opts"/"rates"/"weak"/"bg"/"nucl" to match the Python-side tags
 * (e.g. cpr_log(cfg, "bg", "Background a(t,T) ready in %.2f s", dt) prints
 * "[bg-c] Background a(t,T) ready in 0.42 s"). */
void cpr_log(const CPRConfig *cfg, const char *tag, const char *fmt, ...);

/* True when stdout can carry the decorative UTF-8 output (banner box, rule).
 * A Windows console in a legacy code page -- cp1252, the default outside UTF-8
 * mode -- renders those bytes as mojibake, so callers fall back to ASCII; a
 * redirected stdout is a plain byte stream and takes UTF-8 unchanged. Mirrors
 * main.py's console_encodable(), which guards the same two places on the
 * Python side (where the mismatch raises UnicodeEncodeError rather than
 * garbling). */
int cpr_console_takes_utf8(void);

#endif /* CPRIMAT_LOG_H */
