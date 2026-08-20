/* detailed_balance.csv: csv_split's five-field path. */
#include "fuzz.h"
#include "network_data.h"
#include <stdlib.h>

int LLVMFuzzerTestOneInput(const uint8_t *data, size_t size)
{
    const char *path = fuzz_write_file("db.csv", data, size);
    CPRDetailedBalanceTable t;
    char *err = NULL;
    if (cpr_load_detailed_balance(path, &t, &err) == 0)
        cpr_detailed_balance_free(&t);
    else
        free(err);
    return 0;
}

int LLVMFuzzerInitialize(int *argc, char ***argv)
{
    (void)argc; (void)argv;
    return 0;
}
