## What and why

<!-- What does this change, and why? Link the issue: Closes #123 -->

## Standard / clause

<!-- Standard, edition and clause or equation number the change implements or fixes (n/a for tooling). -->

## Checklist

- [ ] Formulas cite the equation and clause; names are descriptive.
- [ ] Tests compare against a hand calculation, a table of the standard or a limit case.
- [ ] `ruff check .` and `ruff format --check src tests` are clean; `pytest` is green.
- [ ] If `tests/golden/reference.json` changed: every changed number is explained in the commit and in
      `KNOWN_ISSUES.md`.
- [ ] Out-of-range inputs raise `ValueError` instead of returning an unsupported number.
- [ ] `CHANGELOG.md` and `README.md` updated when behaviour or scope changed.
- [ ] No scans, PDFs or long excerpts of copyrighted standards are included.
- [ ] Conventional commit messages.
