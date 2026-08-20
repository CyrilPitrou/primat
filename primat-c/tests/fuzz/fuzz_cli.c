/* Command-line parsing: cpr_cli_main's two argv passes, cpr_parse_literal and
 * cpr_config_set_by_name's did-you-mean formatter.
 *
 * The input is split on NUL into argv words. Every invocation ends with
 * --list-reactions, the last subcommand to return before the solve, so both
 * argv passes and cpr_config_validate run on every input while no run takes
 * minutes. --data_dir is pinned first to a sandbox root whose csv/, nuclear/
 * and NEVO/ are symlinks and whose cache tree is a private empty copy, so a
 * fuzzed --cache-clear cannot touch the shipped caches. */
#include "fuzz.h"
#include "cli.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <unistd.h>

#define MAX_ARGS 64
#define MAX_ARGLEN 256

static char g_root[1024];

int LLVMFuzzerInitialize(int *argc, char ***argv)
{
    (void)argc; (void)argv;
    snprintf(g_root, sizeof(g_root), "%s/root", fuzz_tmpdir());
    char cmd[4096];
    snprintf(cmd, sizeof(cmd),
             "mkdir -p '%s/cache_plasma_weak/weak' '%s/cache_plasma_weak/plasma' && "
             "for d in csv nuclear NEVO; do ln -sfn \"$(cd '%s' && pwd)/$d\" '%s'/$d; done",
             g_root, g_root, fuzz_data_dir(), g_root);
    if (system(cmd) != 0) fprintf(stderr, "fuzz_cli: sandbox setup failed\n");
    /* The CLI is chatty; its output is not what is under test. */
    if (!freopen("/dev/null", "w", stdout)) return 0;
    if (!freopen("/dev/null", "w", stderr)) return 0;
    return 0;
}

int LLVMFuzzerTestOneInput(const uint8_t *data, size_t size)
{
    char *argv[MAX_ARGS];
    static char store[MAX_ARGS][MAX_ARGLEN];
    int argc = 0;

    argv[argc++] = (char *)"primat-c";
    argv[argc++] = (char *)"--data_dir";
    argv[argc++] = g_root;

    size_t pos = 0;
    while (pos < size && argc < MAX_ARGS - 2) {
        size_t end = pos;
        while (end < size && data[end] != '\0') end++;
        size_t n = end - pos;
        if (n > MAX_ARGLEN - 1) n = MAX_ARGLEN - 1;
        memcpy(store[argc], data + pos, n);
        store[argc][n] = '\0';
        argv[argc] = store[argc];
        argc++;
        pos = end + 1;
    }
    argv[argc++] = (char *)"--list-reactions";
    argv[argc] = NULL;

    cpr_cli_main(argc, argv);
    return 0;
}
