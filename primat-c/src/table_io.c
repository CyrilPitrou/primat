/* table_io.c -- see table_io.h. The one numeric-table reader every loader in the
 * port goes through: whitespace- or comma-separated columns, `#` comments,
 * column count either fixed by the caller or detected from the first data row. */

#include "table_io.h"
#include "xalloc.h"

#include <ctype.h>
#include <errno.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Splits `line` in place on whitespace and/or commas, returning pointers to
 * up to `max_fields` field starts in `fields[]` (each NUL-terminated by
 * overwriting the separator). Returns the field count. */
static size_t split_fields(char *line, char **fields, size_t max_fields)
{
    size_t n = 0;
    char *p = line;
    while (*p && n < max_fields) {
        while (*p && (isspace((unsigned char)*p) || *p == ',')) p++;
        if (!*p) break;
        fields[n++] = p;
        while (*p && !isspace((unsigned char)*p) && *p != ',') p++;
        if (*p) { *p = '\0'; p++; }
    }
    return n;
}

int cpr_read_line(FILE *f, char *buf, size_t bufsize)
{
    if (!fgets(buf, (int)bufsize, f))
        return 0;
    size_t n = strlen(buf);
    if (n == 0 || buf[n - 1] == '\n')
        return 1;
    /* No newline: either the line is longer than the buffer, or the file ends
     * without one. Only the first is a truncation. */
    int ch;
    int truncated = 0;
    while ((ch = fgetc(f)) != EOF) {
        truncated = 1;
        if (ch == '\n')
            break;
    }
    return truncated ? -1 : 1;
}

int cpr_table_read(const char *path, size_t n_cols_hint, CPRTable *out,
                    char **errmsg)
{
    FILE *f = fopen(path, "r");
    if (!f) {
        char buf[4352];
        snprintf(buf, sizeof(buf), "cannot open table file '%s'", path);
        *errmsg = strdup(buf);
        return 1;
    }

    out->cols = NULL;
    out->n_cols = n_cols_hint;
    out->n_rows = 0;
    size_t cap = 0;

    char line[8192];
    int lineno = 0;
    char *fields[256];

    int rc;
    while ((rc = cpr_read_line(f, line, sizeof(line))) != 0) {
        lineno++;
        char *s = line;
        while (isspace((unsigned char)*s)) s++;
        if (rc < 0) {
            /* A comment or blank line may be any length -- numpy.loadtxt reads
             * the shipped tables that way and skips it whole, so skipping it
             * here keeps the two backends agreeing. A truncated *data* line
             * cannot be parsed, and must not be parsed in part. */
            if (*s == '\0' || *s == '#')
                continue;
            char buf[4352];
            snprintf(buf, sizeof(buf),
                      "%s:%d: data line longer than %zu characters",
                      path, lineno, sizeof(line) - 1);
            *errmsg = strdup(buf);
            fclose(f);
            cpr_table_free(out);
            return 1;
        }
        if (*s == '\0' || *s == '#')
            continue;

        char linecopy[8192];
        strncpy(linecopy, line, sizeof(linecopy) - 1);
        linecopy[sizeof(linecopy) - 1] = '\0';
        size_t nf = split_fields(linecopy, fields, 256);
        if (nf == 0)
            continue;

        if (out->n_cols == 0) {
            out->n_cols = nf;
        } else if (nf != out->n_cols) {
            char buf[4352];
            snprintf(buf, sizeof(buf),
                      "%s:%d: expected %zu columns, found %zu",
                      path, lineno, out->n_cols, nf);
            *errmsg = strdup(buf);
            fclose(f);
            cpr_table_free(out);
            return 1;
        }

        if (out->n_rows == cap) {
            cap = cap ? cap * 2 : 64;
            if (!out->cols) {
                out->cols = CPR_XCALLOC(out->n_cols, sizeof(double *));
                for (size_t c = 0; c < out->n_cols; c++)
                    out->cols[c] = NULL;
            }
            for (size_t c = 0; c < out->n_cols; c++)
                out->cols[c] = CPR_XREALLOC(out->cols[c], cap * sizeof(double));
        }

        for (size_t c = 0; c < out->n_cols; c++) {
            errno = 0;
            char *endptr;
            double v = strtod(fields[c], &endptr);
            /* ERANGE alone is not necessarily an error: strtod also sets
             * it on underflow to a subnormal (e.g. a rate-table column
             * value like 1.617129e-308, seen in
             * Li7_d__Li8_p_primat.txt) -- a successfully parsed, merely
             * tiny, value. Only overflow to +-HUGE_VAL is a genuine
             * parse failure. */
            if (endptr == fields[c] || (errno == ERANGE && fabs(v) == HUGE_VAL)) {
                char buf[4352];
                snprintf(buf, sizeof(buf),
                          "%s:%d: cannot parse '%s' as a number",
                          path, lineno, fields[c]);
                *errmsg = strdup(buf);
                fclose(f);
                cpr_table_free(out);
                return 1;
            }
            out->cols[c][out->n_rows] = v;
        }
        out->n_rows++;
    }

    fclose(f);

    /* No data rows is an error whatever the caller's column hint was. Checking
     * n_cols too let an empty (or comments-only) file through whenever a hint
     * was given, and every caller then indexes cols[c][0]. */
    if (out->n_rows == 0) {
        char buf[4352];
        snprintf(buf, sizeof(buf), "%s: no data rows found", path);
        *errmsg = strdup(buf);
        cpr_table_free(out);
        return 1;
    }

    return 0;
}

void cpr_table_free(CPRTable *t)
{
    if (!t->cols) return;
    for (size_t c = 0; c < t->n_cols; c++)
        free(t->cols[c]);
    free(t->cols);
    t->cols = NULL;
    t->n_cols = t->n_rows = 0;
}
