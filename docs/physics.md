# Physics

primat solves two coupled problems. The **background** fixes how fast the
universe expands and cools: the photon and per-flavour neutrino temperatures
through non-instantaneous decoupling, QED corrections to the plasma equation
of state, and the resulting `a(t)`, `T(t)` and Hubble rate. The **nuclear
network** then integrates the reactions that build the light elements against
that history, with the n↔p weak rates — Born plus radiative, finite-mass,
thermal and spectral-distortion corrections — setting the neutron fraction the
whole outcome hinges on.

Two documents describe that physics in full. Neither is reproduced here.

**The Physics Reports paper** is the primary reference for every formula in
the code, and what to cite:

> Pitrou, Coc, Uzan, Vangioni, *Physics Reports* **04** (2018) 005
> ([arXiv:1801.08023](https://arxiv.org/abs/1801.08023)).

Source comments cite it by equation number; an annotated copy lives at
`biblio/Pitrou_etal_PhysReptArxivVersion.pdf` in the repository.

**The LaTeX manual** is this package's own write-up — usage, plasma
thermodynamics, weak interactions, nuclear reactions, sensitivity, and
appendices A–G including the full reaction and nuclide tables. It is built
from the current code (its figures are regenerated from the public API, its
reaction tables from the shipped network data), so it is the place to look for
what *this implementation* does rather than what the method paper describes.
It is a PDF in the repository, under
[`manual/`](https://github.com/CyrilPitrou/primat/tree/master/manual), and
`manual/README.md` explains how to rebuild it.

For the code-level view of the same material, {doc}`extending` names the three
extension points and the modules behind them, and the {doc}`api/index` pages
carry each module's own physics notes.
