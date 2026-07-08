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
