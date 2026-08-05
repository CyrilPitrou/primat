#include "ini.h"

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static char *trim(char *s)
{
    while (isspace((unsigned char)*s)) s++;
    if (*s == '\0') return s;
    char *end = s + strlen(s) - 1;
    while (end > s && isspace((unsigned char)*end)) *end-- = '\0';
    return s;
}

int cpr_ini_load(CPRConfig *cfg, const char *path, CPRParamList *collect,
                 char **errmsg)
{
    FILE *f = fopen(path, "r");
    if (!f) {
        char buf[4352];
        snprintf(buf, sizeof(buf), "cannot open ini file '%s'", path);
        *errmsg = strdup(buf);
        return 1;
    }

    char line[4096];
    int lineno = 0;
    while (fgets(line, sizeof(line), f)) {
        lineno++;
        char *s = trim(line);
        if (*s == '\0' || *s == '#' || *s == ';')
            continue;

        /* Split at the first '=' if present, else the first run of
         * whitespace (the "KEY VALUE" form). */
        char *eq = strchr(s, '=');
        char *key, *val;
        if (eq) {
            *eq = '\0';
            key = trim(s);
            val = trim(eq + 1);
        } else {
            char *sp = s;
            while (*sp && !isspace((unsigned char)*sp)) sp++;
            if (*sp == '\0') {
                fprintf(stderr, "%s:%d: warning: ignoring line with no value: '%s'\n",
                        path, lineno, s);
                continue;
            }
            *sp = '\0';
            key = trim(s);
            val = trim(sp + 1);
        }
        if (*key == '\0') {
            fprintf(stderr, "%s:%d: warning: ignoring line with empty key\n", path, lineno);
            continue;
        }

        /* An empty value is not a literal (cpr_parse_literal's closing
         * comment): "network =" must not silently become the empty string. */
        if (*val == '\0') {
            char buf[CPR_PARAM_VAL_LEN];
            snprintf(buf, sizeof(buf), "%s:%d: key '%.200s' has an empty value",
                     path, lineno, key);
            *errmsg = strdup(buf);
            fclose(f);
            return 1;
        }

        CPRParam value = cpr_parse_literal(val);
        char *set_err = NULL;
        int rc = cpr_config_set_by_name(cfg, key, value, &set_err);
        if (rc == CPR_SET_OK) {
            /* Retain it for the MC workers (see ini.h); cpr_paramlist_add
             * copies both halves, so `value` pointing into cpr_parse_literal's
             * static scratch is fine. */
            if (collect)
                cpr_paramlist_add(collect, key, value);
        } else if (rc == CPR_SET_UNKNOWN_KEY && !cfg->strict_params) {
            /* Warn and ignore -- PRIMATConfig's strict_params=False default. */
            fprintf(stderr, "%s:%d: warning: %s\n", path, lineno,
                    set_err ? set_err : "could not set key");
        } else {
            /* A type mismatch (always), or an unknown key under
             * strict_params=True: both raise on the Python side, so both are
             * fatal here. Continuing past a type mismatch is also what used to
             * leave a freed string in the config. */
            char buf[CPR_PARAM_VAL_LEN];
            snprintf(buf, sizeof(buf), "%s:%d: %s%s", path, lineno,
                     set_err ? set_err : "could not set key",
                     rc == CPR_SET_UNKNOWN_KEY ? " [strict_params=True]" : "");
            *errmsg = strdup(buf);
            free(set_err);
            fclose(f);
            return 1;
        }
        free(set_err);
    }

    fclose(f);
    return 0;
}
