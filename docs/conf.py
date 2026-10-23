"""Sphinx configuration of the python-gearbox documentation.

Build it with ``sphinx-build -W -b html docs docs/_build/html`` from the root of the repository
(``pip install -e ".[docs]"`` first).
"""

import os
import sys
from datetime import date

sys.path.insert(0, os.path.abspath("../src"))

import gearbox  # noqa: E402

project = "python-gearbox"
author = "Eduardo M. Fírvida Donéstevez"
copyright = f"2015-{date.today().year}, {author}"
version = release = gearbox.__version__

extensions = [
    "myst_parser",
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.napoleon",
    "sphinx.ext.intersphinx",
    "sphinx.ext.mathjax",
    "sphinx.ext.viewcode",
    "sphinx_copybutton",
    "sphinx_design",
    "sphinx_sitemap",
    "sphinxext.opengraph",
]

source_suffix = {".md": "markdown", ".rst": "restructuredtext"}
root_doc = "index"
templates_path = ["_templates"]
# the standards (PDF and scanned pages) and the work folders of the figures are not part of the docs
exclude_patterns = ["_build", "scribd", "figures/work", "figures/cache", "**/__pycache__"]

# -- MyST ------------------------------------------------------------------------------------------
myst_enable_extensions = ["amsmath", "attrs_inline", "colon_fence", "deflist", "dollarmath"]
myst_heading_anchors = 3
myst_dmath_double_inline = True

# -- API reference ---------------------------------------------------------------------------------
autosummary_generate = True
autosummary_imported_members = False
autodoc_default_options = {
    "members": True,
    "show-inheritance": True,
    "member-order": "bysource",
}
autodoc_typehints = "description"
autodoc_typehints_description_target = "documented_params"
autodoc_class_signature = "separated"
autodoc_preserve_defaults = True
autodoc_mock_imports = ["gmsh", "matplotlib"]
napoleon_google_docstring = True
napoleon_numpy_docstring = False
napoleon_use_rtype = False
napoleon_attr_annotations = True
add_module_names = False
python_use_unqualified_type_names = True
# the docstrings use ``x`` for symbols and code, as in the README
default_role = "code"

intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "numpy": ("https://numpy.org/doc/stable/", None),
    "matplotlib": ("https://matplotlib.org/stable/", None),
}

# -- HTML ------------------------------------------------------------------------------------------
html_theme = "pydata_sphinx_theme"
html_title = f"python-gearbox {release}"
html_static_path = ["_static"]
html_css_files = ["custom.css"]
html_logo = "_static/logo.svg"
html_favicon = "_static/logo.svg"
html_show_sourcelink = False
html_theme_options = {
    "github_url": "https://github.com/efirvida/python-gearbox",
    "icon_links": [
        {
            "name": "Zenodo DOI",
            "url": "https://doi.org/10.5281/zenodo.3660526",
            "icon": "fa-solid fa-quote-right",
        },
    ],
    "logo": {"text": "python-gearbox"},
    "navbar_align": "left",
    "header_links_before_dropdown": 5,
    "show_toc_level": 2,
    "navigation_with_keys": False,
    "use_edit_page_button": True,
    "footer_start": ["copyright"],
    "footer_end": ["sphinx-version", "theme-version"],
}
html_context = {
    "github_user": "efirvida",
    "github_repo": "python-gearbox",
    "github_version": "master",
    "doc_path": "docs",
}

# -- Search engines and social cards ------------------------------------------------------------
html_baseurl = "https://efirvida.github.io/python-gearbox/"
html_extra_path = ["robots.txt"]
sitemap_url_scheme = "{link}"
sitemap_excludes = ["search.html", "genindex.html"]
ogp_site_url = html_baseurl
ogp_site_name = "python-gearbox"
ogp_image = "https://raw.githubusercontent.com/efirvida/python-gearbox/master/docs/figures/mesh_pair_3d.png"
ogp_description_length = 200
ogp_enable_meta_description = True
