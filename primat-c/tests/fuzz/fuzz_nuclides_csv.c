/* nuclides.csv: config.c's load_nuclides, reached the only way it can be --
 * through cpr_config_init_defaults on a data root holding the fuzzed file. */
#include "fuzz.h"
#include "config.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static char g_root[1024];

int LLVMFuzzerInitialize(int *argc, char ***argv)
{
    (void)argc; (void)argv;
    /* The root needs nuclear/ next to csv/ for cpr_config_validate; only
     * csv/nuclides.csv is fuzzed. */
    fuzz_write_file("root/nuclear/.keep", (const uint8_t *)"", 0);
    snprintf(g_root, sizeof(g_root), "%s/root", fuzz_tmpdir());
    return 0;
}

int LLVMFuzzerTestOneInput(const uint8_t *data, size_t size)
{
    fuzz_write_file("root/csv/nuclides.csv", data, size);
    CPRConfig cfg;
    char *err = NULL;
    if (cpr_config_init_defaults(&cfg, g_root, &err) != 0)
        free(err);
    cpr_config_free(&cfg);
    return 0;
}
