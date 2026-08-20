/* Cache files: the fingerprint-header reader plus the three-column data read
 * every cache consumer performs on a header match (weak_rates.c:1288). */
#include "fuzz.h"
#include "cache.h"
#include "table_io.h"
#include <stdlib.h>

int LLVMFuzzerTestOneInput(const uint8_t *data, size_t size)
{
    const char *path = fuzz_write_file("cache.txt", data, size);
    char *hash = cpr_cache_read_fingerprint_hash(path);
    free(hash);
    CPRTable t;
    char *err = NULL;
    if (cpr_table_read(path, 3, &t, &err) == 0)
        cpr_table_free(&t);
    else
        free(err);
    return 0;
}

int LLVMFuzzerInitialize(int *argc, char ***argv)
{
    (void)argc; (void)argv;
    return 0;
}
