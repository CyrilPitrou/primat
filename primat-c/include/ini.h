/* ini.h -- minimal ".ini"-ish config file loader (port of cli.py's --set
 * semantics' .ini handling).
 *
 * Format: one "KEY=VALUE" or "KEY VALUE" per (trimmed) line; lines starting
 * with '#' or ';', and blank lines, are ignored. No sections. VALUE is
 * parsed exactly like --set (cpr_parse_literal): int, then double, then
 * true/false/none (case-insensitive), else the literal string (surrounding
 * quotes stripped if present).
 */
#ifndef CPRIMAT_INI_H
#define CPRIMAT_INI_H

#include "config.h"

/* Loads `path`, applying every KEY=VALUE line to `cfg` via
 * cpr_config_set_by_name (so p_<rxn>/delta_<rxn> keys and ordinary
 * fields are both handled, with the same type-coercion rules).
 *
 * Error policy, mirroring the Python side exactly (see config.h's
 * CPR_SET_* codes):
 *   - unknown key   -> warning on stderr, line ignored, load continues
 *                      (PRIMATConfig's strict_params=False default). When
 *                      cfg->strict_params is set, it is fatal instead, as
 *                      PRIMATConfig raises for strict_params=True.
 *   - wrong type    -> fatal (PRIMATConfig raises TypeError unconditionally).
 *   - unreadable file -> fatal.
 *
 * `collect` (may be NULL) receives a copy of every pair successfully applied,
 * in file order. The CLI passes its Monte-Carlo override list here: MC workers
 * rebuild their config from defaults + that list, so anything the ini set but
 * the list omits silently vanishes from every sample (the sigmas would then
 * describe a different model than the central value printed beside them).
 *
 * Returns 0 on success, nonzero with *errmsg set (caller frees) otherwise. */
int cpr_ini_load(CPRConfig *cfg, const char *path, CPRParamList *collect,
                 char **errmsg);

#endif /* CPRIMAT_INI_H */
