# wheels/ — do not delete

This directory holds a single committed wheel,
`primat-<version>-cp312-*-linux_x86_64.whl`, that is **load-bearing for the
public demo at [primat.streamlit.app](https://primat.streamlit.app)**.

## The deployment chain

1. Streamlit Community Cloud installs the demo app's dependencies from the
   repo-root `requirements.txt`.
2. That file's last line points at `./wheels/primat-<version>-cp312-*-linux_x86_64.whl`
   (a relative path into this directory) instead of `pip install primat`,
   so the demo always runs the exact commit's code rather than whatever is
   currently on PyPI.
3. The wheel itself is produced by `.github/workflows/build_linux.yml`
   (`workflow_dispatch`, manually triggered) and committed here by hand.

None of these three pieces (this wheel, the `requirements.txt` line, or
`build_linux.yml`) is legacy or redundant with the PyPI `wheels.yml`
publishing workflow — deleting any of them breaks the website. See
`CLAUDE.md`'s "Streamlit Cloud deployment chain" section for the full
picture.

## Keeping this in sync with a version bump

Whenever `pyproject.toml`'s `version` is bumped:

1. Run `build_linux.yml` (`workflow_dispatch`) to produce the new wheel.
2. Commit the new wheel file under `wheels/`.
3. Update the filename referenced in `requirements.txt` to match.
4. Delete the old wheel file (a stale one left behind just wastes repo space
   -- `requirements.txt` only ever points at one at a time).

`tests/test_docs_consistency.py` pins that the wheel filename quoted in
`requirements.txt` (a) exists in this directory and (b) has the same version
as `pyproject.toml`, so a missed step above fails a test instead of silently
leaving the website on an old build.
