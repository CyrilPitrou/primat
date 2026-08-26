/* cli.h -- the `primat-c` executable's entry point.
 *
 * cpr_cli_main() is main() in all but name: it parses argv (and any --ini
 * file) into a validated CPRConfig, then either services a standalone
 * subcommand that needs no solve -- --help, --version, --credits,
 * --list-params, --list-reactions, --cache-info, --cache-clear -- or runs the
 * full BBN computation via cprimat_run (api.h), optionally with Monte-Carlo
 * uncertainty propagation (mc.h), and prints and writes the results.
 *
 * Returns the process exit status: 0 on success, 2 on any failure -- a usage
 * error, a rejected configuration, or a failed solve. Nothing here returns 1.
 */
#ifndef CPRIMAT_CLI_H
#define CPRIMAT_CLI_H

int cpr_cli_main(int argc, char **argv);

#endif /* CPRIMAT_CLI_H */
