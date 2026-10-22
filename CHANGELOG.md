# Changelog

All notable changes are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and the project uses
[Semantic Versioning](https://semver.org/). Numerical changes are explained in more detail in
[KNOWN_ISSUES.md](KNOWN_ISSUES.md).

## [Unreleased]

## [0.2.0] - unreleased

### Added
* Double helical gears (`Gear.double_helical`) in ISO 6336 (KI-38).
* Internal gears (`Gear.internal`) in ISO 6336-1/2/3 and in the AGMA 2001/908 pitting calculation (KI-39).
* ISO 6336-1 mesh load factor, methods A and C of the dynamic factor, static strength, limited pitting,
  grinding notches (`Y_Sg`), shaper-cutter root geometry; AGMA 2001-D04 allowable stresses and safety
  factors; ISO 1328-1:2013 flank tolerances.
* Contact ratio limited by the base circle, ISO 6336-2 formula (33) (KI-40).

### Changed
* Complete refactor: Python >= 3.11, immutable typed dataclasses, descriptive names, test suite.
* ISO 6336 updated from the 2006 to the **2019** edition. **Numbers differ from 0.1.x**; see
  [KNOWN_ISSUES.md](KNOWN_ISSUES.md) (KI-01 to KI-37) for every change and its effect.
* AGMA geometry factors `I` and `J` verified against the AGMA 908-B89 tables.

### Not implemented
See "Limitations and pending work" in the [README](README.md).

## [1.0.0]

Original release (ISO 6336:2006 and AGMA 2101, spur and helical gears).
