# Your first BBN run

A guided walk from a fresh install to a result you can defend: what primat
computes, what each number it prints means, and how to change one thing and
see the effect. No prior experience with a BBN code assumed. The ten worked
examples in {doc}`index` pick up where this leaves off.

If you have not installed primat yet, {doc}`../installation` is one line.

## 1. Run the standard case

```bash
primat
```

That solves Big Bang Nucleosynthesis at primat's default settings — the
Planck 2018 + BAO baryon density, three Standard-Model neutrino flavours,
every correction on — and prints:

```text
────────────────────────────────────────────────────
          PRIMAT results at T = 0.001 MeV
────────────────────────────────────────────────────
Neff       = 3.04397730
YP (BBN)   = 0.24699907
YP (CMB)   = 0.24567276
He4/H      = 8.2011454e-02
D/H        = 2.4358767e-05
He3/H      = 1.0399348e-05
He3/He4    = 1.2680361e-04
Li7/H      = 5.557664e-10
```

`T = 0.001 MeV` is where the network stops: about 20 minutes in, nuclear
reactions have frozen out and these abundances no longer change.

## 2. Read the numbers

| Line | What it is |
|---|---|
| `Neff` | the relativistic energy density beyond photons, in units of one instantaneously-decoupled neutrino flavour. 3.044 rather than 3 because decoupling is not instantaneous — that is a *prediction* here, not an input |
| `YP (BBN)` | the primordial helium-4 **mass** fraction. The headline BBN number, compared against ~0.245 from extragalactic HII regions |
| `YP (CMB)` | the same helium, in the convention CMB analyses use |
| `D/H` | deuterium per hydrogen, **by number**. The precision cosmology probe: measured to ~1 % in quasar absorption systems |
| `He3/H`, `He3/He4`, `Li7/H` | the remaining light elements. `Li7/H` is the famous lithium problem — predicted ~3× the observed value |

Two conventions to be clear about before you compare anything with a paper:

- **`YP` is a mass fraction; everything else is a number ratio.** `YP` is
  `4 × Y_He4`, where `Y_i = n_i/n_b` is an abundance per baryon.
- **`D/H` is a number ratio**, `n_D/n_H`, which is what observers quote.

The {doc}`../glossary` has a line for each of these, and
{doc}`../howto/output` documents the full result dict the API returns.

## 3. The same run from Python

```python
from primat.backend import run_bbn

result = run_bbn({"Omegabh2": 0.02242})

print(f"YP  (BBN) = {result['YPBBN']:.8f}")  # 0.24699907
print(f"D/H       = {result['DoH']:.7e}")    # 2.4358767e-05
```

`run_bbn` picks the compiled C engine when it is available and the pure-Python
implementation otherwise, and returns the same dict either way (they agree to
about 7e-06 relative on D/H — {doc}`../api/backend`). It also prints a
one-line progress indicator, `[primat]  HT.  MT.  LT.  done.`, naming the
three temperature eras it integrates across.

## 4. Change one thing

The baryon density is the parameter BBN constrains best, so vary it:

```python
for obh2 in (0.021, 0.02242, 0.023):
    r = run_bbn({"Omegabh2": obh2})
    print(f"Omegabh2={obh2}:  YP={r['YPBBN']:.8f}  D/H={r['DoH']:.7e}")
```

```text
Omegabh2=0.021:    YP=0.24636495  D/H=2.7104310e-05
Omegabh2=0.02242:  YP=0.24699907  D/H=2.4358767e-05
Omegabh2=0.023:    YP=0.24724450  D/H=2.3355864e-05
```

D/H falls steeply as baryon density rises — more baryons burn more deuterium —
while `YP` barely moves. That contrast is why D/H measures Ω_b h² and helium
measures the expansion rate. {doc}`../howto/sensitivity` turns this into a
table of ∂ln(observable)/∂ln(parameter) for every input at once, and
{doc}`PosteriorBaryons` into a posterior on Ω_b h².

Every parameter can be set the same way, from the API, the CLI
(`primat --Omegabh2 0.021`, or `--set KEY=VALUE` for anything without a
flag) or an INI file. {doc}`../parameters` lists them all.

## 5. Look inside the run

The abundances above are the endpoint of an evolution. To get the whole thing:

```python
result = run_bbn({"output_time_evolution": True})
ev = result["evolution"]

print(len(ev.t), ev.t[0], ev.t[-1])   # 500 rows, 4.6e-04 s to 1.3e+06 s
print(ev.T_gamma[0], ev.T_gamma[-1])  # 40.0 MeV down to 0.001 MeV
print(list(ev.Y))                     # n, p, H2, H3, He3, He4, Li7, Be7
print(ev.Y["He4"][-1] * 4)            # 0.24699907 -- YP again
```

`ev` is an `EvolutionResult`: `t` [s], `a`, `T_gamma` [MeV], `T_nu["e"]` and
its two siblings, and `Y[<nuclide>]`, each a 500-point array. It is in the
result dict whenever the flag is set, *and* written to
`results/output_tables.tsv` — pass `output_file=None` to keep it in memory
only. {doc}`../howto/output` documents the file's schema, and
{doc}`AbundanceEvolution` plots it.

## Where to go next

| If you want to | Go to |
|---|---|
| an uncertainty on these numbers | {doc}`../howto/rate-variation-mc` |
| a bigger reaction network | {doc}`../howto/networks` |
| a non-standard expansion history | {doc}`../howto/backgrounds` |
| to click rather than type | {doc}`../howto/gui` |
| the physics behind all of it | {doc}`../physics` |
