# Installation

python-gearbox needs Python 3.11 or newer. The strength calculations depend only on NumPy, SciPy
and Jinja2.

```bash
pip install "python-gearbox @ git+https://github.com/efirvida/python-gearbox"       # strength calculations
pip install "python-gearbox[fem] @ git+https://github.com/efirvida/python-gearbox"  # plus meshing and plots
```

```{admonition} Version 0.2 is not on PyPI yet
:class: warning
`pip install python-gearbox` still installs the 0.1 release of 2015, which has a different API.
Once 0.2 is published, `pip install python-gearbox` and `pip install "python-gearbox[fem]"` will do.
```

## Finite-element analyses

The `fem` extra installs [Gmsh](https://gmsh.info) (meshing) and matplotlib (plots). Solving needs
**CalculiX**, a separate program:

```bash
conda install -c conda-forge calculix    # or: apt install calculix-ccx
```

Make sure `ccx` is on your `PATH`, or point the `GEARBOX_CCX` environment variable at it.
{py:func}`~gearbox.fem.calculix.find_ccx` searches an explicit path, then `GEARBOX_CCX`, then the
`PATH`. On Linux, Gmsh needs a few system libraries (`libglu1-mesa`, `libxrender1`, `libxcursor1`,
`libxft2`, `libxinerama1`).

## Development setup

```bash
git clone https://github.com/efirvida/python-gearbox
cd python-gearbox
pip install -e ".[dev]"
pytest
```

To build this documentation:

```bash
pip install -e ".[docs]"
sphinx-build -W -b html docs docs/_build/html
```
