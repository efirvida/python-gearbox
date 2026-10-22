"""Rewrite ``tests/golden/export`` from the current code.

Run from the repository root: ``python tests/golden/regenerate_export.py``. Like
``reference.json``, the files must be regenerated in the same commit as any change that
intentionally alters the exported geometry, and the diff explained.
"""

from __future__ import annotations

import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from tests.scenarios import scenario  # noqa: E402
from tests.test_export_golden import DIMENSIONS, GOLDEN, SCENARIOS, TARGETS, _run  # noqa: E402


def main() -> None:
    for name in SCENARIOS:
        transmission = scenario(name)
        for kind in ("pinion", "wheel", "pair"):
            for target in TARGETS:
                for dimension in DIMENSIONS:
                    with tempfile.TemporaryDirectory() as tmp:
                        out = Path(tmp) / "out"
                        written = _run(kind, transmission, target, dimension, out)
                        text = written.read_text().replace(str(out), "<OUTPUT_DIR>")
                        (GOLDEN / name / written.name).write_text(text)


if __name__ == "__main__":
    main()
