# Command-line interface

The `primat` console script exposes the most-used run options directly. The
reference below is generated automatically from the argument parser, so it never
drifts from `primat --help`.

```{eval-rst}
.. argparse::
   :module: primat.cli
   :func: _build_parser
   :prog: primat
```

:::{tip}
Anything not exposed as a named flag can still be set via the repeatable
`--set KEY=VALUE` escape hatch for any `PRIMATConfig` key, e.g.
`primat --set T_end_MeV=1e-4 --set network=large`.
:::

:::{tip}
Per-reaction rate variations use the `p_<reaction>` (log-normal, in units of
the reaction's tabulated 1σ factor) and `delta_<reaction>` (additive) keys.
They depend on the network, so they have no fixed list and no named flags —
`primat --list-reactions` prints the names the selected `--network`/`--amax`
accepts:

```bash
primat --list-reactions --network large --amax 8
primat --set p_d_p__He3_g=1.0        # raise that rate by one sigma
```
:::
