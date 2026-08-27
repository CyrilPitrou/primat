# PyPiGuide.md — Publishing a new `primat` release to PyPI

This is the recurring checklist for shipping a new version, with every
irreversible action flagged and a way to test up to (but not past) each one.

Current repo state: the name is claimed and Trusted Publishing is already
wired up on both indexes. **PyPI already has `0.3.0` and `0.3.1`; TestPyPI
already has `0.3.2`.** There is no first-upload/name-claiming step left to
do — every release from now on follows the same recurring path below.

## Legend

- 🟢 **Reversible** — redo it, undo it, no lasting effect outside your
  own machine/repo.
- 🟡 **Hard to undo** — affects a shared system (GitHub) but can be
  cleaned up; mistakes are recoverable with effort.
- 🔴 **Irreversible** — affects PyPI's permanent, public, append-only
  index. Cannot be undone. Do these last, and only once you've verified
  everything beforehand.

---

## Step 1 — 🟢 Bump the version everywhere

`pyproject.toml` is the single source of truth; every other place repeats it
and none of them updates itself. This table is the whole list — if a file is
not here, it does not carry the version.

| # | File | What to change | Checked by |
|---|---|---|---|
| 1 | `pyproject.toml` | `version` | — (source of truth) |
| 2 | `primat-c/include/config.h` | `CPRIMAT_VERSION` | `test_cprimat_version_matches_pyproject` |
| 3 | `CITATION.cff` | `version` **and** `date-released` (the real date) | `test_citation_cff_version_matches_pyproject` |
| 4 | `manual/` | rename `primat_documentation_vX.Y.Z.tex`/`.pdf`, update the `.tex` title page and intro, and `manual/README.md`'s four references | `test_manual_declares_the_current_version` |
| 5 | `CHANGELOG.md` | turn `[Unreleased]` into `## [X.Y.Z] - <date>` (the real date), and check every user-visible change since the last release is listed | its **date** is cross-checked against row 3; its contents are not machine-checkable |
| 6 | `wheels/` + `requirements.txt` | rebuild and re-point the Streamlit wheel — **Step 3**, not here | `test_streamlit_wheel_matches_pyproject_version` |

Rows 2–4 fail a test if you skip them, so `pytest tests/test_docs_consistency.py`
is the check for this step. Row 6 is Step 3 rather than here because it needs a
CI run. Only row 5's *contents* are left to you.

All of this is a normal commit — reversible.

## Step 2 — 🟢 Build and check locally, no network upload

Before touching any external service, verify the artifacts you'd
eventually publish actually work:

```bash
# sdist + wheel for your current platform only
python -m build

# Validate package metadata (long_description renders, no malformed
# classifiers, etc.) without uploading anywhere
pip install twine
twine check dist/*

# Sanity-install into a throwaway venv from the built wheel, not -e .
python -m venv /tmp/primat-check
source /tmp/primat-check/bin/activate
pip install dist/primat-*.whl
python -c "from primat import PRIMAT; print(PRIMAT({}).solve()['DoH'])"
deactivate
```

You can repeat this as many times as you want. Nothing here touches
PyPI, TestPyPI, or GitHub.

### Optional: exercise the multi-platform wheel matrix locally

`cibuildwheel` can run on your own machine before any of it goes near
GitHub Actions:

```bash
pip install cibuildwheel
cibuildwheel --platform macos   # builds the matrix for the OS you're on
```

This is the closest you can get to rehearsing `wheels.yml`'s
`build_wheels` job with zero network publish step and zero GitHub
involvement. It won't catch a Windows-specific MSVC regression (you're
not on Windows) or the `aarch64` QEMU cross-build, but it does validate
the `setup.py`/`pyproject.toml` packaging metadata and the `primat-c`
extension build flags.

---

## Step 3 — 🟢 Rebuild and commit the Streamlit demo wheel

The public demo at **primat.streamlit.app** installs `primat` from the
committed wheel under `wheels/` (see `wheels/README.md` for the full
deployment chain). Every version bump needs this wheel refreshed:

1. Trigger `.github/workflows/build_linux.yml` (`workflow_dispatch`) to
   produce the new `primat-X.Y.Z-cp312-*-linux_x86_64.whl`.
2. Commit it under `wheels/`, deleting the previous version's wheel.
3. Update the filename in `requirements.txt`'s last line to point at it.

Skipping this step doesn't block the PyPI release below, but it leaves the
Streamlit demo silently serving the old version.

---

## Step 4 — 🟢 Dry-run on TestPyPI

1. Trigger `wheels.yml` via `workflow_dispatch` with `publish_testpypi:
   true` → wheels + sdist get built and uploaded to **test**.pypi.org
   (`publish-testpypi` job, `environment: testpypi`, Trusted Publishing
   already configured — no credentials to enter).
2. Verify end-to-end: `pip install -i https://test.pypi.org/simple/ primat`
   into a clean venv, run the validation script
   (`runfiles/primat_run.py`), and confirm the result matches the
   documented tolerances in `tests/README.md`'s "Validation reference".

Why 🟢 this time (unlike the first-ever upload): the `(name, version)` pair
being uploaded is the new version you're about to release, so there is no
prior TestPyPI upload to collide with — a failed or mistaken attempt just
means bumping to `X.Y.Zrc1`, `rc2`, etc. and retrying.

---

## Step 5 — 🔴 Tag `vX.Y.Z`, publish the GitHub release, let `wheels.yml` upload to real PyPI

The GitHub release's `published` event triggers `wheels.yml`'s full
pipeline against the real `pypi` index, via the dedicated `publish-pypi`
job (`environment: pypi`, Trusted Publishing already registered), which
only ever fires on `release: published` and is unreachable from
`workflow_dispatch`.

Irreversible because:
- Once `primat-X.Y.Z` (wheels + sdist) lands on real PyPI, that exact
  version's files can never be replaced — only deleted (hiding it from
  new installs, but not erasing it from anyone who already resolved it,
  and not freeing the version number for reuse with different content).
- Anyone in the world can `pip install primat==X.Y.Z` from the moment
  the first wheel finishes uploading. There is no "private" undo.

Pre-flight checklist (everything above must already be true):
- [ ] Step 1: every row of its table done — `pytest tests/test_docs_consistency.py`
      green, `CHANGELOG.md` dated and complete.
- [ ] Step 2: local build + `twine check` clean.
- [ ] Step 3: Streamlit demo wheel rebuilt and committed (or explicitly
      deferred — see Step 3's note, not a release blocker).
- [ ] Step 4: full TestPyPI install-and-run dry run matches the published
      tolerances.
- [ ] The version string in `pyproject.toml` is exactly what you intend
      to ship — this is your last chance to change it before it's
      permanent.
- [ ] You (not an agent) create the `vX.Y.Z` git tag and the GitHub
      release, since this is the action that actually fires the
      irreversible publish — treat the "Publish release" button on
      GitHub as the point of no return.

After this, the only thing left to verify is `pip install primat` from
a clean machine/venv with no `-i test.pypi.org` flag, confirming the
real index serves what you expect.

---

## Step 6 — 🟡 Zenodo archival DOI (per-release housekeeping)

`CITATION.cff` (repo root) gives the repo a "Cite this repository" button
on GitHub and structured metadata for GitHub/Zenodo/other citation
tooling. The GitHub↔Zenodo integration is a one-time toggle already done;
what's recurring per release:

1. Keep `CITATION.cff`'s `version`/`date-released` fields in sync with
   `pyproject.toml` at each release (mirroring the `CPRIMAT_VERSION` sync
   habit).
2. The GitHub release from Step 5 automatically triggers Zenodo to
   snapshot that tagged version, mint a version-specific DOI (the
   version-independent "concept DOI" in `README.md`'s badge always
   resolves to it — no README edit needed per release).

🟡 because a Zenodo record, once minted for a published release, is meant
to be permanent (that is the point of an archival DOI), even though the
CITATION.cff edit itself is a normal reversible commit.
