/* decays.txt: cpr_load_decays' sscanf row parser. */
#include "fuzz.h"
#include "network_data.h"
#include <stdlib.h>

int LLVMFuzzerTestOneInput(const uint8_t *data, size_t size)
{
    const char *path = fuzz_write_file("decays.txt", data, size);
    CPRDecayTable t;
    char *err = NULL;
    if (cpr_load_decays(path, &t, &err) == 0)
        cpr_decay_table_free(&t);
    else
        free(err);
    return 0;
}

int LLVMFuzzerInitialize(int *argc, char ***argv)
{
    (void)argc; (void)argv;
    return 0;
}
