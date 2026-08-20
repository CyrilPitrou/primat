/* Rate tables: cpr_table_read's auto-width path plus the validator every
 * caller runs on its columns (network_data.c:1299). */
#include "fuzz.h"
#include "table_io.h"
#include "network_data.h"
#include <stdlib.h>

int LLVMFuzzerTestOneInput(const uint8_t *data, size_t size)
{
    const char *path = fuzz_write_file("rate.txt", data, size);
    CPRTable t;
    char *err = NULL;
    if (cpr_table_read(path, 0, &t, &err) == 0) {
        if (t.n_cols >= 2) {
            const double *e = t.n_cols >= 3 ? t.cols[2] : NULL;
            char *verr = NULL;
            if (cpr_validate_rate_table(t.cols[0], t.cols[1], e, t.n_rows,
                                        "fuzz", &verr) != 0)
                free(verr);
        }
        cpr_table_free(&t);
    } else {
        free(err);
    }
    return 0;
}

int LLVMFuzzerInitialize(int *argc, char ***argv)
{
    (void)argc; (void)argv;
    return 0;
}
