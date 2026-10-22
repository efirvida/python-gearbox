# Golden data

- `baseline_original.json` – numbers of the original 0.1.x code (legacy key names). Historical, never edited, not compared by tests.
- `generate_baseline_original_api.py`, `generate_export_original_api.py` – scripts that produced the original baseline / export files with the old API.
- `reference.json` – current expected numbers (new field names) for every scenario in `tests/scenarios.py`: gear and mesh geometry, ISO 6336 / AGMA 2101 results, profile-shift optimisation. Initialised from the original numbers.
- `regenerate_reference.py` – rebuilds `reference.json` from the current code.
- `export/` – golden files for the export tests.

Regenerate:

    .venv/bin/python tests/golden/regenerate_reference.py

Rule: a numerical fix must regenerate `reference.json` in the same commit, and the commit message must explain every changed number in the diff.

`regenerate_export.py` rewrites `tests/golden/export/` from the current code; use it (and explain the diff) when a change intentionally alters the exported geometry.
