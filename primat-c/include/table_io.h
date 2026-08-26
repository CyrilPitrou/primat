/* table_io.h -- generic reader for the whitespace/comma-separated numeric
 * tables used throughout data/ (per-reaction rate tables, NEVO CSVs, QED
 * tables, electron_thermo_<hash>.txt). All of these share the same shape:
 * zero or more '#'-prefixed header/comment lines, then one row of N
 * numbers per data line, columns separated by whitespace and/or commas
 * (commas appear in the NEVO files, whitespace in the rate/QED tables).
 *
 * This module is pure I/O (no physics, no resampling) -- it hands back a
 * column-major double** plus the row count; per-format semantics (which
 * column is T9, what the header documents, ...) live in network_data.c/
 * neutrino_history.c/qed_pressure.c, which call this and interpret the
 * result.
 */
#ifndef CPRIMAT_TABLE_IO_H
#define CPRIMAT_TABLE_IO_H

#include <stddef.h>
#include <stdio.h>

typedef struct {
    double **cols;   /* cols[c][r], c in [0, n_cols), r in [0, n_rows) */
    size_t n_cols;
    size_t n_rows;
} CPRTable;

/* Reads every non-comment, non-blank line of `path` as `n_cols` numbers
 * (auto-detected from the first data line if n_cols_hint == 0; every
 * subsequent data line must then match that column count exactly, or this
 * fails -- a malformed/truncated row is treated as a hard error, not
 * silently dropped, since these files drive physics results). Returns 0 on
 * success (caller must cpr_table_free the result), nonzero with *errmsg set
 * (caller frees) on failure (missing file, inconsistent column count,
 * unparsable number). */
int cpr_table_read(const char *path, size_t n_cols_hint, CPRTable *out,
                    char **errmsg);

/* Releases the column arrays and zeroes the struct; `t` itself is the
 * caller's. */
void cpr_table_free(CPRTable *t);

/* Reads one whole line from `f` into `buf` (at most `bufsize`-1 bytes plus the
 * NUL), with the trailing newline kept, exactly as fgets does.
 *
 * The difference is what happens to a line that does not fit: fgets returns
 * the first chunk and hands the remainder back on the next call, as a line of
 * its own. Every reader here is line-oriented, so that turned the tail of an
 * over-long comment into a data row -- silently, and only for the one input
 * that straddles the buffer. This consumes the rest of the line instead and
 * says so.
 *
 * Returns 1 on a complete line, 0 at end of file, and -1 when the line was
 * too long (the remainder is consumed, so the caller may report and continue
 * from the next line). */
int cpr_read_line(FILE *f, char *buf, size_t bufsize);

#endif /* CPRIMAT_TABLE_IO_H */
