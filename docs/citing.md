# Citing primat

If you use primat in published work, please cite the Physics Reports paper:

> Pitrou, Coc, Uzan, Vangioni, *Physics Reports* **04** (2018) 005
> ([arXiv:1801.08023](https://arxiv.org/abs/1801.08023)).

```bibtex
@article{Pitrou:2018cgg,
    author        = "Pitrou, Cyril and Coc, Alain and Uzan, Jean-Philippe and Vangioni, Elisabeth",
    title         = "{Precision big bang nucleosynthesis with improved Helium-4 predictions}",
    journal       = "Phys. Rept.",
    volume        = "754",
    pages         = "1--66",
    year          = "2018",
    eprint        = "1801.08023",
    archivePrefix = "arXiv",
    primaryClass  = "astro-ph.CO",
    doi           = "10.1016/j.physrep.2018.04.005"
}
```

The same entry is available at runtime as `primat.__citation__` (see
`primat/credits.py`), so scripts can print it without hard-coding it.

## Citing the software itself

If you want to cite the specific `primat` software/release you used (as
distinct from the physics methods paper above), the repository ships a
[`CITATION.cff`](https://github.com/CyrilPitrou/primat/blob/master/CITATION.cff)
at its root, which GitHub renders as a "Cite this repository" button on the
repo page. Once Zenodo–GitHub archival is enabled (see `PyPiGuide.md` Step
6), each tagged release also gets its own archival DOI through Zenodo.

