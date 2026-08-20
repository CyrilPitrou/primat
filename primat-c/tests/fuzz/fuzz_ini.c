/* The .ini loader: cpr_ini_load + every cpr_config_set_by_name path it can
 * reach, including cpr_parse_literal's type dispatch. */
#include "fuzz.h"
#include "config.h"
#include "ini.h"
#include <stdlib.h>
#include <string.h>

int LLVMFuzzerTestOneInput(const uint8_t *data, size_t size)
{
    const char *path = fuzz_write_file("run.ini", data, size);
    CPRConfig cfg;
    char *err = NULL;
    if (cpr_config_init_defaults(&cfg, fuzz_data_dir(), &err) != 0) {
        free(err);
        cpr_config_free(&cfg);
        return 0;
    }
    CPRParamList cp;
    memset(&cp, 0, sizeof(cp));
    if (cpr_ini_load(&cfg, path, &cp, &err) != 0) {
        free(err);
    } else {
        char *verr = NULL;
        if (cpr_config_validate(&cfg, &verr) != 0) free(verr);
    }
    cpr_paramlist_free(&cp);
    cpr_config_free(&cfg);
    return 0;
}

int LLVMFuzzerInitialize(int *argc, char ***argv)
{
    (void)argc; (void)argv;
    return 0;
}
