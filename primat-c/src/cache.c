#include "cache.h"
#include "constants.h"
#include "xalloc.h"

#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "compat_posix.h"  /* unistd.h / getpid, portable across POSIX & MSVC */

/* ===========================================================================
 * SHA-256 (public-domain style, single-buffer one-shot implementation --
 * no streaming API needed since every fingerprint JSON blob is short).
 * ===========================================================================
 */
static const uint32_t SHA256_K[64] = {
    0x428a2f98,0x71374491,0xb5c0fbcf,0xe9b5dba5,0x3956c25b,0x59f111f1,0x923f82a4,0xab1c5ed5,
    0xd807aa98,0x12835b01,0x243185be,0x550c7dc3,0x72be5d74,0x80deb1fe,0x9bdc06a7,0xc19bf174,
    0xe49b69c1,0xefbe4786,0x0fc19dc6,0x240ca1cc,0x2de92c6f,0x4a7484aa,0x5cb0a9dc,0x76f988da,
    0x983e5152,0xa831c66d,0xb00327c8,0xbf597fc7,0xc6e00bf3,0xd5a79147,0x06ca6351,0x14292967,
    0x27b70a85,0x2e1b2138,0x4d2c6dfc,0x53380d13,0x650a7354,0x766a0abb,0x81c2c92e,0x92722c85,
    0xa2bfe8a1,0xa81a664b,0xc24b8b70,0xc76c51a3,0xd192e819,0xd6990624,0xf40e3585,0x106aa070,
    0x19a4c116,0x1e376c08,0x2748774c,0x34b0bcb5,0x391c0cb3,0x4ed8aa4a,0x5b9cca4f,0x682e6ff3,
    0x748f82ee,0x78a5636f,0x84c87814,0x8cc70208,0x90befffa,0xa4506ceb,0xbef9a3f7,0xc67178f2
};

static uint32_t rotr32(uint32_t x, int n) { return (x >> n) | (x << (32 - n)); }

static void sha256(const unsigned char *msg, size_t len, unsigned char out[32])
{
    uint32_t h[8] = {
        0x6a09e667,0xbb67ae85,0x3c6ef372,0xa54ff53a,
        0x510e527f,0x9b05688c,0x1f83d9ab,0x5be0cd19
    };

    /* Pad: msg || 0x80 || zeros || 64-bit big-endian bit length, to a
     * multiple of 64 bytes. */
    size_t padded_len = ((len + 9 + 63) / 64) * 64;
    unsigned char *buf = CPR_XCALLOC(padded_len, 1);
    memcpy(buf, msg, len);
    buf[len] = 0x80;
    uint64_t bitlen = (uint64_t)len * 8;
    for (int i = 0; i < 8; i++)
        buf[padded_len - 1 - i] = (unsigned char)(bitlen >> (8 * i));

    for (size_t off = 0; off < padded_len; off += 64) {
        uint32_t w[64];
        for (int i = 0; i < 16; i++)
            w[i] = ((uint32_t)buf[off + 4*i] << 24) | ((uint32_t)buf[off + 4*i + 1] << 16)
                 | ((uint32_t)buf[off + 4*i + 2] << 8) | (uint32_t)buf[off + 4*i + 3];
        for (int i = 16; i < 64; i++) {
            uint32_t s0 = rotr32(w[i-15], 7) ^ rotr32(w[i-15], 18) ^ (w[i-15] >> 3);
            uint32_t s1 = rotr32(w[i-2], 17) ^ rotr32(w[i-2], 19) ^ (w[i-2] >> 10);
            w[i] = w[i-16] + s0 + w[i-7] + s1;
        }

        uint32_t a=h[0],b=h[1],c=h[2],d=h[3],e=h[4],f=h[5],g=h[6],hh=h[7];
        for (int i = 0; i < 64; i++) {
            uint32_t S1 = rotr32(e,6) ^ rotr32(e,11) ^ rotr32(e,25);
            uint32_t ch = (e & f) ^ (~e & g);
            uint32_t t1 = hh + S1 + ch + SHA256_K[i] + w[i];
            uint32_t S0 = rotr32(a,2) ^ rotr32(a,13) ^ rotr32(a,22);
            uint32_t maj = (a & b) ^ (a & c) ^ (b & c);
            uint32_t t2 = S0 + maj;
            hh=g; g=f; f=e; e=d+t1; d=c; c=b; b=a; a=t1+t2;
        }
        h[0]+=a; h[1]+=b; h[2]+=c; h[3]+=d; h[4]+=e; h[5]+=f; h[6]+=g; h[7]+=hh;
    }
    free(buf);

    for (int i = 0; i < 8; i++) {
        out[4*i]   = (unsigned char)(h[i] >> 24);
        out[4*i+1] = (unsigned char)(h[i] >> 16);
        out[4*i+2] = (unsigned char)(h[i] >> 8);
        out[4*i+3] = (unsigned char)(h[i]);
    }
}

char *cpr_sha256_hex16(const char *json_str)
{
    unsigned char digest[32];
    sha256((const unsigned char *)json_str, strlen(json_str), digest);
    char *hex = CPR_XMALLOC(17);
    for (int i = 0; i < 8; i++)
        snprintf(hex + 2*i, 3, "%02x", digest[i]);
    return hex;
}

/* ===========================================================================
 * Canonical JSON serialisation -- must match
 * json.dumps(d, sort_keys=True, separators=(",", ":")) exactly.
 * ===========================================================================
 */

/* Python float repr: the shortest decimal string that round-trips to the
 * same IEEE-754 double. We brute-force the precision (1..17 significant
 * digits via "%.*e"/"%.*g") rather than porting Grisu/dtoa -- the
 * fingerprint fields are always "nice" config values (0.0, 0.001, 40.0,
 * ...), so this is never on a hot path and correctness-by-construction
 * (verify round-trip with strtod) matters more than speed here. */
/* Python's float repr (CPython's format_float_short, mode 'r'): find the
 * shortest decimal digit string that round-trips to the same double, then
 * format it fixed-point if its decimal exponent is in [-4, 16), else
 * scientific -- e.g. repr(40.0)=="40.0", repr(0.001)=="0.001",
 * repr(1e+16)=="1e+16", repr(1e-05)=="1e-05". Using "%.*e" (always
 * scientific) to find the shortest round-tripping digit string avoids %g's
 * premature switch to scientific notation at low precision (e.g. "%.1g" of
 * 40.0 is "4e+01", not "40"), which is the bug this replaces. */
static void format_python_float(double v, char *buf, size_t bufsize)
{
    if (v == 0.0) {
        snprintf(buf, bufsize, "%s", signbit(v) ? "-0.0" : "0.0");
        return;
    }

    char e_str[64];
    int p;
    for (p = 1; p <= 17; p++) {
        snprintf(e_str, sizeof(e_str), "%.*e", p - 1, v);
        if (strtod(e_str, NULL) == v)
            break;
    }

    /* e_str is "[-]d[.ddd]e±EE". Extract sign, significant digits (no '.'),
     * and the decimal exponent of the leading digit. */
    const char *s = e_str;
    int neg = 0;
    if (*s == '-') { neg = 1; s++; }
    char digits[32];
    int nd = 0;
    digits[nd++] = *s++; /* leading digit */
    if (*s == '.') {
        s++;
        while (*s && *s != 'e' && *s != 'E')
            digits[nd++] = *s++;
    }
    while (*s && *s != 'e' && *s != 'E') s++;
    s++; /* skip 'e' */
    int decexp = atoi(s);

    char out[64];
    size_t len = 0;
    if (neg) out[len++] = '-';

    if (decexp >= -4 && decexp < 16) {
        if (decexp >= 0) {
            int int_digits = decexp + 1;
            for (int i = 0; i < int_digits; i++)
                out[len++] = (i < nd) ? digits[i] : '0';
            out[len++] = '.';
            if (int_digits < nd) {
                for (int i = int_digits; i < nd; i++) out[len++] = digits[i];
            } else {
                out[len++] = '0';
            }
        } else {
            out[len++] = '0';
            out[len++] = '.';
            for (int i = 0; i < -decexp - 1; i++) out[len++] = '0';
            for (int i = 0; i < nd; i++) out[len++] = digits[i];
        }
    } else {
        out[len++] = digits[0];
        if (nd > 1) {
            out[len++] = '.';
            for (int i = 1; i < nd; i++) out[len++] = digits[i];
        }
        len += snprintf(out + len, sizeof(out) - len, "e%+03d", decexp);
    }
    out[len] = '\0';
    snprintf(buf, bufsize, "%s", out);
}

static void json_escape_append(char **buf, size_t *cap, size_t *len, const char *s)
{
    size_t need = *len + strlen(s) * 2 + 4;
    if (need > *cap) {
        *cap = need * 2;
        *buf = realloc(*buf, *cap);
    }
    (*buf)[(*len)++] = '"';
    for (const char *p = s; *p; p++) {
        unsigned char c = (unsigned char)*p;
        if (c == '"' || c == '\\') {
            (*buf)[(*len)++] = '\\';
            (*buf)[(*len)++] = (char)c;
        } else if (c == '\n') { (*buf)[(*len)++] = '\\'; (*buf)[(*len)++] = 'n'; }
        else if (c == '\t') { (*buf)[(*len)++] = '\\'; (*buf)[(*len)++] = 't'; }
        else if (c < 0x20) {
            *len += snprintf(*buf + *len, 8, "\\u%04x", c);
        } else {
            (*buf)[(*len)++] = (char)c;
        }
    }
    (*buf)[(*len)++] = '"';
}

static void buf_append(char **buf, size_t *cap, size_t *len, const char *s)
{
    size_t slen = strlen(s);
    if (*len + slen + 1 > *cap) {
        *cap = (*len + slen + 1) * 2;
        *buf = realloc(*buf, *cap);
    }
    memcpy(*buf + *len, s, slen + 1);
    *len += slen;
}

char *cpr_fingerprint_json(const CPRFPField *fields, size_t n)
{
    /* Sort a local copy of the field indices by key (byte-wise, matching
     * Python's default string comparison for plain ASCII identifiers). */
    size_t *order = CPR_XMALLOC(n * sizeof(size_t));
    for (size_t i = 0; i < n; i++) order[i] = i;
    for (size_t i = 1; i < n; i++) {
        size_t j = i;
        while (j > 0 && strcmp(fields[order[j-1]].key, fields[order[j]].key) > 0) {
            size_t t = order[j-1]; order[j-1] = order[j]; order[j] = t;
            j--;
        }
    }

    size_t cap = 256, len = 0;
    char *buf = CPR_XMALLOC(cap);
    buf_append(&buf, &cap, &len, "{");
    for (size_t k = 0; k < n; k++) {
        const CPRFPField *f = &fields[order[k]];
        if (k > 0) buf_append(&buf, &cap, &len, ",");
        json_escape_append(&buf, &cap, &len, f->key);
        buf_append(&buf, &cap, &len, ":");
        char num[64];
        switch (f->value.type) {
        case CPR_NONE:
            buf_append(&buf, &cap, &len, "null");
            break;
        case CPR_BOOL:
            buf_append(&buf, &cap, &len, f->value.v.b ? "true" : "false");
            break;
        case CPR_INT:
            snprintf(num, sizeof(num), "%ld", f->value.v.i);
            buf_append(&buf, &cap, &len, num);
            break;
        case CPR_DOUBLE:
            format_python_float(f->value.v.d, num, sizeof(num));
            buf_append(&buf, &cap, &len, num);
            break;
        case CPR_STRING:
            json_escape_append(&buf, &cap, &len, f->value.v.s);
            break;
        }
    }
    buf_append(&buf, &cap, &len, "}");
    free(order);
    return buf;
}

char *cpr_fingerprint_hash(const CPRFPField *fields, size_t n)
{
    char *json = cpr_fingerprint_json(fields, n);
    char *hash = cpr_sha256_hex16(json);
    free(json);
    return hash;
}

/* ===========================================================================
 * Weak-rate / thermal fingerprint builders (port of weak_rates/cache.py).
 * ===========================================================================
 */
static CPRParam pb(int b) { CPRParam p; p.type = CPR_BOOL; p.v.b = b; return p; }
static CPRParam pi(long i) { CPRParam p; p.type = CPR_INT; p.v.i = i; return p; }
static CPRParam pd(double d) { CPRParam p; p.type = CPR_DOUBLE; p.v.d = d; return p; }
static CPRParam ps(const char *s) {
    CPRParam p;
    if (s) { p.type = CPR_STRING; p.v.s = s; }
    else { p.type = CPR_NONE; }
    return p;
}

/* ---------------------------------------------------------------------------
 * Physical-constants hash -- the C mirror of cache_utils.constants_hash().
 *
 * Hashes every field of g_const (26 of them, exactly the fields of Python's
 * `Constants` dataclass, which is what dataclasses.asdict() serialises), so
 * that any edit to a physical constant invalidates every cache computed with
 * the old value. See cache.h for the contract and cache_utils.constants_hash's
 * docstring for the design rationale.
 *
 * The @property-derived quantities (sW2, MeV_to_Kelvin, mB, ...) are
 * deliberately absent on both sides: they are pure functions of the fields,
 * so they carry no extra information, and keeping them out means the hash
 * never depends on a float the C compiler might evaluate differently from
 * CPython (e.g. cpr_erg() divides by `second` where the Python property
 * divides by `second**2` -- inert at second = 1, but it would poison a hash
 * that included it).
 * ------------------------------------------------------------------------- */
const char *cpr_constants_hash(void)
{
    /* Memoised for the process: the constants are frozen, so this is the C
     * equivalent of Python's @lru_cache(maxsize=1). 17 = 16 hex digits + NUL. */
    static char cached[17];
    if (cached[0] != '\0')
        return cached;

    /* Idempotent; makes the helper safe to call before any explicit init. */
    cpr_constants_init();

    /* Listed in the Python dataclass's declaration order for side-by-side
     * comparison with constants.py; cpr_fingerprint_json sorts by key anyway,
     * exactly as json.dumps(sort_keys=True) does. */
    const CPRFPField fields[] = {
        /* ---- CGS base units (natural units: all 1) ---- */
        {"Kelvin",    pd(g_const.Kelvin)},
        {"second",    pd(g_const.second)},
        {"cm",        pd(g_const.cm)},
        {"gram",      pd(g_const.gram)},
        /* ---- Fundamental constants ---- */
        {"kB",        pd(g_const.kB)},
        {"clight",    pd(g_const.clight)},
        {"hbar",      pd(g_const.hbar)},
        {"Mpc",       pd(g_const.Mpc)},
        {"MeV",       pd(g_const.MeV)},
        {"keV",       pd(g_const.keV)},
        /* ---- Electroweak sector ---- */
        {"alphaem",   pd(g_const.alphaem)},
        {"GF",        pd(g_const.GF)},
        {"mZ",        pd(g_const.mZ)},
        /* ---- Fermion masses [MeV] ---- */
        {"me",        pd(g_const.me)},
        {"mn",        pd(g_const.mn)},
        {"mp",        pd(g_const.mp)},
        /* ---- CMB ---- */
        {"T0CMB",     pd(g_const.T0CMB)},
        /* ---- Weak-rate nuclear-structure constants ---- */
        {"gA",        pd(g_const.gA)},
        {"kappa_p",   pd(g_const.kappa_p)},
        {"kappa_n",   pd(g_const.kappa_n)},
        {"Vud",       pd(g_const.Vud)},
        {"radproton", pd(g_const.radproton)},
        /* ---- Atomic masses ---- */
        {"ma",        pd(g_const.ma)},
        {"He4Overma", pd(g_const.He4Overma)},
        {"HOverma",   pd(g_const.HOverma)},
        /* ---- Standard-model effective neutrino number ---- */
        {"Neff_SM",   pd(g_const.Neff_SM)},
    };
    /* Compile-time guard: if a field is added to CPRConstants (and to Python's
     * Constants) but not to the array above, the two hashes silently diverge
     * and every shared cache file becomes a cross-backend miss. Keep this
     * count in step with constants.py's field list. */
    _Static_assert(sizeof(fields) / sizeof(fields[0]) == 26,
                   "cpr_constants_hash must hash all 26 Constants fields "
                   "(see primat/constants.py); update both on any change");

    char *hash = cpr_fingerprint_hash(fields, sizeof(fields) / sizeof(fields[0]));
    snprintf(cached, sizeof(cached), "%s", hash);
    free(hash);
    return cached;
}

/* WEAK_RATE_FORMAT_VERSION in weak_rates/cache.py (see its changelog comment
 * for what each generation changed; v4 pays off the v2/v3 bumps that were
 * documented but never applied, and adds the munuOverTnu / nevo_grid_file /
 * sampling_temperature_per_decade fields below; v5 adds constants_hash to both
 * fingerprints, so that editing any physical constant invalidates the caches
 * computed with the old value instead of silently reloading them). */
#define WEAK_RATE_FORMAT_VERSION 5

size_t cpr_weak_rate_fingerprint(const CPRConfig *cfg, CPRFPField *out)
{
    size_t n = 0;
    out[n++] = (CPRFPField){"format_version", pi(WEAK_RATE_FORMAT_VERSION)};
    out[n++] = (CPRFPField){"sampling_nTOp_per_decade", pi(cfg->sampling_nTOp_per_decade)};
    out[n++] = (CPRFPField){"radiative_corrections", pb(cfg->radiative_corrections)};
    out[n++] = (CPRFPField){"finite_mass_corrections", pb(cfg->finite_mass_corrections)};
    /* Physical constants (v5): the rate integrands read me, alphaem, mn, mp,
     * gA, Vud, radproton, kappa_n/p and GF directly. Mirrors Python
     * _weak_rate_fingerprint's "constants_hash" entry; the pointer is to
     * cpr_constants_hash's process-lifetime static buffer, so it safely
     * outlives this call. */
    out[n++] = (CPRFPField){"constants_hash", ps(cpr_constants_hash())};
    /* Effective ξ_e under the historical "munuOverTnu" key: only nu_e
     * shifts the n<->p rates, and using the effective xi_e (= munuOverTnu when
     * munuOverTnu_e is unset) keeps the default-run hash unchanged so shipped
     * data/weak/ caches stay valid. Mirrors Python _weak_rate_fingerprint. */
    out[n++] = (CPRFPField){"munuOverTnu", pd(cpr_config_xi_nu_e(cfg))};
    out[n++] = (CPRFPField){"QED_corrections", pb(cfg->QED_corrections)};
    out[n++] = (CPRFPField){"incomplete_decoupling", pb(cfg->incomplete_decoupling)};
    out[n++] = (CPRFPField){"spectral_distortions", pb(cfg->spectral_distortions)};
    out[n++] = (CPRFPField){"analytic_distortions", pb(cfg->analytic_distortions)};
    out[n++] = (CPRFPField){"y_SZ", pd(cfg->y_SZ)};
    out[n++] = (CPRFPField){"y_gray", pd(cfg->y_gray)};
    out[n++] = (CPRFPField){"T_start_cosmo_MeV", pd(cfg->T_start_cosmo_MeV)};
    out[n++] = (CPRFPField){"T_end_MeV", pd(cfg->T_end_MeV)};
    /* Sets the node spacing of the linear T_nu(T_gamma) interpolant every rate
     * integrand reads, hence the rates themselves (~1e-5 at the default 600
     * points/decade, ~1e-3 at 40) -- see the matching note in
     * weak_rates/cache.py's _WEAK_RATE_BG_FIELDS. */
    out[n++] = (CPRFPField){"sampling_temperature_per_decade",
                            pi(cfg->sampling_temperature_per_decade)};
    out[n++] = (CPRFPField){"nevo_file", ps(cfg->nevo_file)};
    out[n++] = (CPRFPField){"nevo_spectral_file", ps(cfg->nevo_spectral_file)};
    /* Pairs with nevo_spectral_file: grid nodes + distortion columns jointly
     * define the dFDneu the spectral-distortion term integrates. */
    out[n++] = (CPRFPField){"nevo_grid_file", ps(cfg->nevo_grid_file)};
    out[n++] = (CPRFPField){"nevo_file_prefix", ps(cfg->nevo_file_prefix)};
    /* custom_background mode takes the Tg grid behind the T_nu(T_gamma)
     * interpolant from the table file's OWN T range (cpr_bg_init_custom), not
     * from T_start_cosmo_MeV/T_end_MeV -- so nothing above distinguishes one
     * custom table from another, and two different backgrounds shared a single
     * cached nTOp table. Keyed on the path, as nevo_file is. Emitted only when
     * set, so a run without custom_background hashes byte-identically to
     * before and keeps hitting the shipped caches (same conditional trick as
     * munuOverTnu). Mirrors weak_rates/cache.py's _weak_rate_fingerprint. */
    if (cfg->custom_background && cfg->custom_background[0] != '\0')
        out[n++] = (CPRFPField){"custom_background", ps(cfg->custom_background)};
    return n; /* 19 entries, +1 iff custom_background is set;
                 sampling_nTOp_per_decade/radiative_corrections/
                 finite_mass_corrections each appear once here already
                 (the Python dict's apparent "duplicate" assignment from
                 looping over _WEAK_RATE_BG_FIELDS after the literal dict
                 is a no-op re-write of the same key -- not represented
                 twice in a dict, so not duplicated here either). */
}

/* T_end_MeV is deliberately NOT in this fingerprint: L_CCRTh_compute clamps
 * the integral to exactly 0 below ~10^8.2 K regardless of cfg.T_end, and the
 * thermal cache grid is now built down to that fixed floor (CCRTH_T_MIN in
 * weak_rates.c) rather than down to cfg.T_end (mirrors
 * weak_rates/cache.py's _THERMAL_BG_FIELDS). Including it here only forced
 * spurious cache misses -- and the multi-minute recompute that goes with
 * them -- whenever a run changed T_end_MeV alone. */
size_t cpr_thermal_fingerprint(const CPRConfig *cfg, CPRFPField *out)
{
    size_t n = 0;
    out[n++] = (CPRFPField){"format_version", pi(WEAK_RATE_FORMAT_VERSION)};
    out[n++] = (CPRFPField){"sampling_nTOp_thermal_per_decade", pi(cfg->sampling_nTOp_thermal_per_decade)};
    /* Physical constants (v5): the CCRTh correction is itself O(alphaem) and
     * its integrands carry me. Mirrors Python _thermal_fingerprint. */
    out[n++] = (CPRFPField){"constants_hash", ps(cpr_constants_hash())};
    out[n++] = (CPRFPField){"T_start_cosmo_MeV", pd(cfg->T_start_cosmo_MeV)};
    out[n++] = (CPRFPField){"QED_corrections", pb(cfg->QED_corrections)};
    out[n++] = (CPRFPField){"incomplete_decoupling", pb(cfg->incomplete_decoupling)};
    out[n++] = (CPRFPField){"nevo_file", ps(cfg->nevo_file)};
    out[n++] = (CPRFPField){"nevo_file_prefix", ps(cfg->nevo_file_prefix)};
    /* Effective ξ_e, same historical key as the weak-rate fingerprint: the
     * thermal integrands' Chitilde carries exp(znu*(en - sgnq*q) - sgnq*xi_nu)
     * (th_chitilde in weak_rates.c), so the CCRTh table is xi_e-specific --
     * measured ~4e-3 of the base rate at xi_e = 0.3, T = 1e10 K. Mirrors
     * weak_rates/cache.py's _thermal_fingerprint. */
    out[n++] = (CPRFPField){"munuOverTnu", pd(cpr_config_xi_nu_e(cfg))};
    return n; /* 9 entries */
}

/* ===========================================================================
 * Cache file read/write (port of cache_utils.read_cache_fingerprint_hash /
 * write_cache_with_fingerprint).
 * ===========================================================================
 */
char *cpr_cache_read_fingerprint_hash(const char *path)
{
    FILE *f = fopen(path, "r");
    if (!f) return NULL;
    char line[2048];
    char *result = NULL;
    while (fgets(line, sizeof(line), f)) {
        if (line[0] != '#') break;
        const char *prefix = "# fingerprint_hash:";
        if (strncmp(line, prefix, strlen(prefix)) == 0) {
            char *val = line + strlen(prefix);
            while (*val == ' ') val++;
            char *end = val + strlen(val);
            while (end > val && (end[-1] == '\n' || end[-1] == '\r' || end[-1] == ' '))
                end--;
            *end = '\0';
            result = strdup(val);
            break;
        }
    }
    fclose(f);
    return result;
}

int cpr_cache_write(const char *path, const CPRFPField *fields, size_t n_fields,
                     const char *col_header, double **columns, size_t n_cols,
                     size_t n_rows, const char *provenance)
{
    char *json = cpr_fingerprint_json(fields, n_fields);
    char *hash = cpr_sha256_hex16(json);

    char tmp_path[4200];
    snprintf(tmp_path, sizeof(tmp_path), "%s.tmp.%d", path, (int)getpid());

    FILE *f = fopen(tmp_path, "w");
    if (!f) { free(json); free(hash); return 1; }

    if (col_header && *col_header)
        fprintf(f, "# %s\n", col_header);
    fprintf(f, "# fingerprint_hash: %s\n", hash);
    fprintf(f, "# fingerprint: %s\n", json);
    if (provenance && *provenance)
        fprintf(f, "# provenance: %s\n", provenance);
    for (size_t r = 0; r < n_rows; r++) {
        for (size_t c = 0; c < n_cols; c++) {
            if (c > 0) fputc(' ', f);
            fprintf(f, "%.18e", columns[c][r]);
        }
        fputc('\n', f);
    }
    fclose(f);
    free(json);
    free(hash);

    if (rename(tmp_path, path) != 0)
        return 1;
    return 0;
}
