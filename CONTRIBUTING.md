# Contributing

Thanks for helping. This library implements published standards, so the bar is **traceability**: every
formula must point at its clause and be checked against a hand calculation or a table of the standard.

## Before you start

* **Bugs and questions:** open an [issue](../../issues/new/choose). For a wrong number, include the
  input, the result, what you expected, and the clause or table that says so.
* **Features:** open an issue first, so the scope (which standard, which clause) is agreed before coding.
* **Security:** see [SECURITY.md](SECURITY.md).
* Read [KNOWN_ISSUES.md](KNOWN_ISSUES.md): it lists what is verified, what is a known deviation and
  what is out of scope.

## Development setup

```bash
git clone https://github.com/efirvida/python-gearbox.git
cd python-gearbox
python -m venv .venv && . .venv/bin/activate
pip install -e ".[dev]"
pytest -q
ruff check . && ruff format --check src tests
```

Python 3.11 or newer. CalculiX is optional (its test is skipped when `ccx` is not installed).

## Rules for changes

1. **Cite the standard.** Name the equation and clause in the docstring or comment
   (`# ISO 6336-2:2019 formula (17)`). Use full descriptive names for variables; the standard's symbol
   goes in the docstring.
2. **Test behaviour, not code.** Compare against a value calculated by hand from the formula, a table of
   the standard, or a limit case (e.g. an internal gear with many teeth tends to the rack case).
3. **Do not move existing numbers silently.** `tests/golden/reference.json` pins the results of the
   reference scenarios. If a change moves it, regenerate it with
   `python tests/golden/regenerate_reference.py`, review the diff, and explain every changed number in
   the commit and in [KNOWN_ISSUES.md](KNOWN_ISSUES.md).
4. **Fail loudly outside the validated range.** Raise `ValueError` with a clear message instead of
   returning a number the standard does not support.
5. **The standards are not distributed.** Do not commit scans, PDFs or long excerpts of copyrighted
   standards; reference them by number and clause.
6. Keep `ruff check` and `ruff format` clean and the test suite green.

## Commits and pull requests

* [Conventional Commits](https://www.conventionalcommits.org/): `feat(iso6336): ...`, `fix(agma2101): ...`,
  `docs: ...`, `test: ...`, `refactor: ...`.
* One logical change per commit, with its tests and docs. Keep pull requests small and focused.
* Fill in the pull request template; link the issue it closes.
* Update [CHANGELOG.md](CHANGELOG.md) under *Unreleased* for user-visible changes.

## Licence

By contributing you agree that your contribution is licensed under the project's
[MIT](LICENSE) licence.
