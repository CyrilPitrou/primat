/* NEVO spectral table + y-grid: a width mismatch between the two is an
 * out-of-bounds read if the loader trusts either width. One fuzzed byte string drives both files, split at the
 * first "%%%" marker. */
#include "fuzz.h"
#include "config.h"
#include "plasma.h"
#include "neutrino_history.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static CPRConfig g_cfg;
static CPRPlasma g_plasma;
static int g_ok;

int LLVMFuzzerInitialize(int *argc, char ***argv)
{
    (void)argc; (void)argv;
    char *err = NULL;
    g_ok = (cpr_config_init_defaults(&g_cfg, fuzz_data_dir(), &err) == 0);
    free(err); err = NULL;
    if (g_ok) {
        g_cfg.spectral_distortions = 1;
        g_cfg.analytic_distortions = 0;
        g_ok = (cpr_plasma_init(&g_plasma, &g_cfg, &err) == 0);
        free(err);
    }
    if (!g_ok) fprintf(stderr, "fuzz_nevo_spectral: setup failed\n");
    return 0;
}

int LLVMFuzzerTestOneInput(const uint8_t *data, size_t size)
{
    if (!g_ok) return 0;
    const uint8_t *sep = NULL;
    for (size_t i = 0; i + 3 <= size; i++)
        if (memcmp(data + i, "%%%", 3) == 0) { sep = data + i; break; }
    size_t n1 = sep ? (size_t)(sep - data) : size;
    const uint8_t *d2 = sep ? sep + 3 : data;
    size_t n2 = sep ? size - n1 - 3 : 0;

    char full[1024], grid[1024];
    snprintf(full, sizeof(full), "%s", fuzz_write_file("nevo_full.csv", data, n1));
    snprintf(grid, sizeof(grid), "%s", fuzz_write_file("nevo_grid.csv", d2, n2));

    free(g_cfg.nevo_spectral_file); g_cfg.nevo_spectral_file = strdup(full);
    free(g_cfg.nevo_grid_file);     g_cfg.nevo_grid_file = strdup(grid);
    CPRNeutrinoHistory nh;
    char *err = NULL;
    if (cpr_neutrino_history_init(&nh, &g_cfg, &g_plasma, &err) == 0)
        cpr_neutrino_history_free(&nh);
    else
        free(err);
    return 0;
}
