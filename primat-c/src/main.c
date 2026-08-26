/* main.c -- the standalone `primat-c` executable, and nothing else. All the
 * work is in cli.c; keeping main() alone here is what lets the unit tests and
 * the Python extension link every other object file without a competing
 * entry point (see the Makefile's filter-out of main.o).
 */
#include "cli.h"

int main(int argc, char **argv)
{
    return cpr_cli_main(argc, argv);
}
