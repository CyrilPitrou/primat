/* test_network_data.c -- loads real rates/nuclear/{networks,data,tables}
 * files and checks counts plus a few known rows, including the small.txt
 * 12-reaction / large.txt 428-reaction counts. */
#include "network_data.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int failures = 0;

#define CHECK(cond, msg) do { \
        if (!(cond)) { printf("FAIL: %s\n", msg); failures++; } \
        else printf("ok: %s\n", msg); \
    } while (0)

static int close(double a, double b) { return fabs(a - b) < 1e-6 * fabs(b) + 1e-300; }

int main(void)
{
    char *err = NULL;

    CPRNetworkList small;
    if (cpr_load_network_list("../primat/data/nuclear/networks/small.txt", &small, &err)) {
        printf("FAIL small.txt: %s\n", err);
        return 1;
    }
    CHECK(small.n == 12, "small.txt has 12 reactions");
    CHECK(strcmp(small.entries[0].name, "n_p__d_g") == 0, "small.txt first reaction name");
    CHECK(strcmp(small.entries[0].table_file, "n_p__d_g_primat.txt") == 0,
          "small.txt first table file");
    cpr_network_list_free(&small);

    CPRNetworkList large;
    if (cpr_load_network_list("../primat/data/nuclear/networks/large.txt", &large, &err)) {
        printf("FAIL large.txt: %s\n", err);
        return 1;
    }
    CHECK(large.n == 428, "large.txt has 428 reactions");
    cpr_network_list_free(&large);

    CPRDecayTable decays;
    if (cpr_load_decays("../primat/data/nuclear/tables/decays.txt", &decays, &err)) {
        printf("FAIL decays.txt: %s\n", err);
        return 1;
    }
    CHECK(decays.n > 0, "decays.txt has rows");
    int found_be7 = 0;
    for (size_t i = 0; i < decays.n; i++) {
        if (strcmp(decays.entries[i].name, "Be7__Li7_Bp") == 0) {
            found_be7 = 1;
            CHECK(close(decays.entries[i].halflife_s, 4.604256e+06), "Be7__Li7_Bp halflife_s");
            CHECK(strstr(decays.entries[i].ref, "Audi") != NULL, "Be7__Li7_Bp ref text");
        }
    }
    CHECK(found_be7, "found Be7__Li7_Bp decay row");
    cpr_decay_table_free(&decays);

    CPRDetailedBalanceTable db;
    if (cpr_load_detailed_balance("../primat/data/csv/detailed_balance.csv", &db, &err)) {
        printf("FAIL detailed_balance.csv: %s\n", err);
        return 1;
    }
    CHECK(db.n > 0, "detailed_balance.csv has rows");
    CHECK(strcmp(db.entries[0].reaction, "B10_He3__C11_d") == 0, "detailed_balance first reaction");
    CHECK(close(db.entries[0].Q_keV, 3196.71), "detailed_balance first Q_keV");
    cpr_detailed_balance_free(&db);

    CPRReactionTable rx;
    if (cpr_load_reactions_large("../primat/data/csv/reactions_large.csv", &rx, &err)) {
        printf("FAIL reactions_large.csv: %s\n", err);
        return 1;
    }
    CHECK(rx.n == 428, "reactions_large.csv has 428 rows");
    CHECK(strcmp(rx.entries[0].reactants, "B10+He3") == 0, "reactions_large first reactants");
    cpr_reaction_table_free(&rx);

    /* Rate-table domain rules, in lockstep with network_data.py's
     * validate_rate_table: an unsorted T9 column used to run to completion and
     * report YP = 0.00000711 against a true 0.24699907. */
    {
        double good_T9[4]   = {1.0e-3, 1.0e-2, 1.0e-1, 1.0};
        double good_rate[4] = {1.0, 2.0, 3.0, 4.0};
        double good_err[4]  = {1.0, 1.0, 1.0, 1.0};
        err = NULL;
        CHECK(cpr_validate_rate_table(good_T9, good_rate, good_err, 4,
                                        "good", &err) == 0,
               "a valid rate table passes");

        struct { const char *label; double T9[4]; double rate[4]; const char *want; }
        bad[] = {
            {"unsorted T9",   {1.0, 1.0e-3, 1.0e-1, 1.0}, {1.0, 2.0, 3.0, 4.0}, "goes backwards"},
            {"duplicate T9",  {1.0e-3, 1.0e-3, 1.0e-1, 1.0}, {1.0, 2.0, 3.0, 4.0}, "repeats"},
            {"zero T9",       {0.0, 1.0e-2, 1.0e-1, 1.0}, {1.0, 2.0, 3.0, 4.0}, "strictly positive"},
            {"negative rate", {1.0e-3, 1.0e-2, 1.0e-1, 1.0}, {1.0, -2.0, 3.0, 4.0}, "must not be negative"},
            {"NaN rate",      {1.0e-3, 1.0e-2, 1.0e-1, 1.0}, {1.0, NAN, 3.0, 4.0}, "non-finite"},
        };
        for (size_t i = 0; i < sizeof(bad) / sizeof(bad[0]); i++) {
            err = NULL;
            int rc = cpr_validate_rate_table(bad[i].T9, bad[i].rate, NULL, 4,
                                               bad[i].label, &err);
            CHECK(rc != 0 && err && strstr(err, bad[i].want) != NULL, bad[i].label);
            free(err);
            err = NULL;
        }

        double one_T9 = 1.0, one_rate = 1.0;
        err = NULL;
        CHECK(cpr_validate_rate_table(&one_T9, &one_rate, NULL, 1, "one row", &err) != 0,
               "a single-row table is refused");
        free(err);
        err = NULL;
    }

    if (failures) {
        printf("%d failure(s)\n", failures);
        return 1;
    }
    printf("all tests passed\n");
    return 0;
}
