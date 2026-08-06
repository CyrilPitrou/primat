/* test_nuclear_network.c -- end-to-end small/large+amax=8 BBN solves
 * (cpr_nuclear_network_solve, the HT->MT->LT era port of
 * primat/nuclear_network.py's NuclearNetwork.solve) against the reference
 * numbers published in tests/README.md's "Validation reference", plus
 * baryon-number conservation at the final state.
 *
 * Tolerance note (read before changing any number below): `rtol` here is far
 * looser than the published +-1e-5 (YP) / +-3e-9 (D/H) bounds, deliberately.
 * BBN abundances are exponentially sensitive to T(t) near weak freeze-out, so
 * the accuracy of the *background* ODEs, not of this network solve, sets how
 * closely the end-to-end numbers land -- which is why background.c decouples
 * their tolerance from cfg->numerical_precision (see BG_ODE_RTOL's comment
 * there). The margin left here still catches every failure this test exists
 * for -- wrong stoichiometry, a wrong reverse-rate cap, a sign error -- by
 * orders of magnitude, while not turning a legitimate small accuracy change
 * into a red suite. The tight bounds are enforced on the Python side, against
 * the same published table. */
#include "constants.h"
#include "plasma.h"
#include "background.h"
#include "network_data.h"
#include "nuclear_network.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int failures = 0;

#define CHECK(cond, msg) do { \
        if (!(cond)) { printf("FAIL: %s\n", msg); failures++; } \
        else printf("ok: %s\n", msg); \
    } while (0)

static int close_rel(double a, double b, double rtol)
{
    return fabs(a - b) <= rtol * fabs(b) + 1e-300;
}

/* Runs one full HT->MT->LT solve for `network`/`amax` (amax<0 => no
 * filter) at the project default config, and checks: it succeeds, baryon
 * number is conserved at the final state to ~1e-9 (Phys. Rep.'s
 * sum_s A_s Y_s = 1 invariant -- a structural check independent of any
 * solver-tolerance question above), and YP(BBN)=4*Y_He4 / D/H / per-
 * nuclide Y match the published reference numbers to within `rtol` (see
 * this file's header comment for why that is looser than the published
 * bounds). `name` is only used in CHECK() messages. */
static void run_and_check(const char *name, const char *network, int amax,
                           double YP_target, double DoH_target,
                           double Yn_target, double Yp_target, double YHe4_target,
                           double rtol)
{
    char *err = NULL;
    CPRConfig cfg;
    if (cpr_config_init_defaults(&cfg, "../primat/data", &err)) {
        printf("FAIL %s: config init: %s\n", name, err); failures++; return;
    }
    free((void *)cfg.network);
    cfg.network = strdup(network);
    if (amax > 0) cfg.amax = amax;

    CPRPlasma pl;
    if (cpr_plasma_init(&pl, &cfg, &err)) {
        printf("FAIL %s: plasma init: %s\n", name, err); failures++; return;
    }
    CPRBackground bg;
    if (cpr_bg_init_standard(&bg, &cfg, &pl, &err)) {
        printf("FAIL %s: background init: %s\n", name, err); failures++; return;
    }
    CPRNuclearRates nr;
    if (cpr_nuclear_rates_init(&nr, &cfg, NULL, &err)) {
        printf("FAIL %s: nuclear rates init: %s\n", name, err); failures++; return;
    }

    CPRNuclearNetwork nn;
    int rc = cpr_nuclear_network_solve(&nn, &cfg, &nr, &bg, &err);
    char msg[160];
    snprintf(msg, sizeof(msg), "%s: cpr_nuclear_network_solve succeeds", name);
    CHECK(rc == 0, msg);
    if (rc) {
        printf("  error: %s\n", err);
        /* Free everything owned so far so LeakSanitizer stays clean even on
         * a solve failure (cfg owns all the init_defaults strdup'd string
         * fields + nuclides.items + the network override strdup'd above). */
        cpr_nuclear_rates_free(&nr);
        cpr_background_free(&bg);
        cpr_plasma_free(&pl);
        cpr_config_free(&cfg);
        return;
    }

    /* Baryon number conservation: sum_s A_s Y_s = 1 (Phys. Rep.) -- the
     * network-correctness invariant. */
    double baryon_sum = 0.0;
    for (size_t i = 0; i < nn.n_species; i++) {
        for (size_t j = 0; j < cfg.nuclides.n; j++)
            if (strcmp(cfg.nuclides.items[j].name, nn.abundance_names[i]) == 0) {
                baryon_sum += (cfg.nuclides.items[j].N + cfg.nuclides.items[j].Z)
                              * nn.Y_final[i];
                break;
            }
    }
    snprintf(msg, sizeof(msg), "%s: baryon number conserved (sum A_s Y_s = 1)", name);
    CHECK(fabs(baryon_sum - 1.0) < 1e-9, msg);

    double Yn = cpr_nuclear_network_get(&nn, "n");
    double Yp = cpr_nuclear_network_get(&nn, "p");
    double YH2 = cpr_nuclear_network_get(&nn, "H2");
    double YHe4 = cpr_nuclear_network_get(&nn, "He4");
    double YPBBN = 4.0 * YHe4;
    double DoH = YH2 / Yp;

    snprintf(msg, sizeof(msg), "%s: YP(BBN) matches the reference within %.0f%%", name, rtol * 100.0);
    CHECK(close_rel(YPBBN, YP_target, rtol), msg);
    snprintf(msg, sizeof(msg), "%s: D/H matches the reference within %.0f%%", name, rtol * 100.0);
    CHECK(close_rel(DoH, DoH_target, rtol), msg);
    snprintf(msg, sizeof(msg), "%s: Yn matches the reference within %.0f%%", name, rtol * 100.0);
    CHECK(close_rel(Yn, Yn_target, rtol), msg);
    snprintf(msg, sizeof(msg), "%s: Yp matches the reference within %.0f%%", name, rtol * 100.0);
    CHECK(close_rel(Yp, Yp_target, rtol), msg);
    snprintf(msg, sizeof(msg), "%s: YHe4 matches the reference within %.0f%%", name, rtol * 100.0);
    CHECK(close_rel(YHe4, YHe4_target, rtol), msg);

    cpr_nuclear_network_free(&nn);
    cpr_nuclear_rates_free(&nr);
    cpr_background_free(&bg);
    cpr_plasma_free(&pl);
    /* cfg owns all the strdup'd default string fields (cpr_config_init_defaults)
     * plus nuclides.items and the network-name override strdup'd above; without
     * this LeakSanitizer (Linux/gcc ASan) reports them all as leaked. */
    cpr_config_free(&cfg);
}

int main(void)
{
    cpr_constants_init();

    /* tests/README.md "Validation reference" numbers (small
     * network, default config: Omegabh2=0.02242, spectral_distortions=
     * QED_corrections=True, numerical_precision=1e-7). 0.1% comfortably
     * covers this port's worst-observed deviation (Yn, the most sensitive
     * trace quantity, ~0.012%) with ~10x margin -- see this file's header
     * comment for how that number was reached. */
    run_and_check("small", "small", -1,
                  0.24700028, 2.43500e-5,
                  3.995347e-16, 7.529409e-01, 6.174973e-02,
                  1.0e-3);

    /* The "large, amax=8" table (the old "medium" network's exact
     * 68-reaction equivalent). */
    run_and_check("large,amax=8", "large", 8,
                  0.24700363, 2.43571e-5,
                  3.994404e-16, 7.529375e-01, 6.175059e-02,
                  1.0e-3);

    if (failures) {
        printf("%d failure(s)\n", failures);
        return 1;
    }
    printf("all tests passed\n");
    return 0;
}
