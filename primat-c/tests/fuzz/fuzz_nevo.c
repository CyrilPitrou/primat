/* NEVO thermo table: cpr_table_read plus build_nevo_table's consumer, reached
 * through cpr_neutrino_history_init with nevo_file pointed at the input. */
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
        /* The spectral branch has its own target; keep this one on the thermo
         * table so its edges are not diluted. */
        g_cfg.spectral_distortions = 0;
        g_ok = (cpr_plasma_init(&g_plasma, &g_cfg, &err) == 0);
        free(err);
    }
    if (!g_ok) fprintf(stderr, "fuzz_nevo: setup failed\n");
    return 0;
}

int LLVMFuzzerTestOneInput(const uint8_t *data, size_t size)
{
    if (!g_ok) return 0;
    const char *path = fuzz_write_file("nevo.csv", data, size);
    free(g_cfg.nevo_file);
    g_cfg.nevo_file = strdup(path);
    CPRNeutrinoHistory nh;
    char *err = NULL;
    if (cpr_neutrino_history_init(&nh, &g_cfg, &g_plasma, &err) == 0)
        cpr_neutrino_history_free(&nh);
    else
        free(err);
    return 0;
}
