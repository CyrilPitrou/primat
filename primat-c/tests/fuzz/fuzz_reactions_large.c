/* reactions_large.csv: the same csv_split path with five string fields. */
#include "fuzz.h"
#include "network_data.h"
#include <stdlib.h>

int LLVMFuzzerTestOneInput(const uint8_t *data, size_t size)
{
    const char *path = fuzz_write_file("rx.csv", data, size);
    CPRReactionTable t;
    char *err = NULL;
    if (cpr_load_reactions_large(path, &t, &err) == 0)
        cpr_reaction_table_free(&t);
    else
        free(err);
    return 0;
}

int LLVMFuzzerInitialize(int *argc, char ***argv)
{
    (void)argc; (void)argv;
    return 0;
}
