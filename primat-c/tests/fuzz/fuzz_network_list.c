/* Network list files under data/nuclear/networks: cpr_load_network_list. */
#include "fuzz.h"
#include "network_data.h"
#include <stdlib.h>

int LLVMFuzzerTestOneInput(const uint8_t *data, size_t size)
{
    const char *path = fuzz_write_file("net.txt", data, size);
    CPRNetworkList list;
    char *err = NULL;
    if (cpr_load_network_list(path, &list, &err) == 0)
        cpr_network_list_free(&list);
    else
        free(err);
    return 0;
}

int LLVMFuzzerInitialize(int *argc, char ***argv)
{
    (void)argc; (void)argv;
    return 0;
}
