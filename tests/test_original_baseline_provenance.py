"""Provenance of the golden reference.

``tests/golden/reference.json`` was initialised from the numbers of the ORIGINAL 0.1.x code
(``tests/golden/baseline_original.json``, legacy key names) after verifying that the refactored
code reproduced them to rel=1e-9. From then on, numerical bug fixes deliberately diverge from the
original, so the baseline is kept only as historical evidence and is NOT compared numerically.
"""

import json
from pathlib import Path

BASELINE = Path(__file__).parent / "golden" / "baseline_original.json"


def test_original_baseline_is_kept_and_valid_json():
    assert BASELINE.is_file()
    assert isinstance(json.loads(BASELINE.read_text()), dict)
