# Configuration file for the Sphinx documentation builder.

from pathlib import Path


def version_from_file(f):
    with f.open("r") as fp:
        return fp.read().strip()


# -- Project information

project = "VAMPyR"
copyright = "TO BE FILLED"
author = "TO BE FILLED"

github_user = "MRChemSoft"
github_repo_name = "vampyr"  # auto-detected from dirname if blank
github_version = "master"

version_file = Path(__file__).absolute().parents[1] / "VERSION"
release = version_from_file(version_file)
version = version_from_file(version_file)

# -- General configuration

extensions = [
    "sphinx.ext.duration",
    "sphinx.ext.doctest",
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.intersphinx",
    "sphinx.ext.napoleon",
    "sphinxcontrib.bibtex",
    "sphinx_thebe",
    "myst_nb",
]

source_suffix = {".rst": "restructuredtext", ".ipynb": "myst-nb", ".myst": "myst-nb"}

# configure sphinxcontrib.bibtex
bibtex_bibfiles = ["bibliography.bib"]

# Settings for myst_nb:
# https://myst-nb.readthedocs.io/en/latest/computation/execute.html
# *only* execute if at least one output is missing.
nb_execution_mode = "off"
nb_execution_show_tb = True

# settings for myst_parser
myst_enable_extensions = [
    "amsmath",
    "colon_fence",
    "deflist",
    "dollarmath",
    "html_image",
]
myst_url_schemes = ("http", "https", "mailto")

autosummary_generate = True

autodoc_default_options = {
    "members": None,
    "imported-members": True,
    "undoc-members": True,
    "show-inheritance": True,
}

intersphinx_mapping = {
    "python": ("https://docs.python.org/3/", None),
    "numpy": ("https://numpy.org/doc/stable/", None),
    "h5py": ("https://docs.h5py.org/en/latest/", None),
    "matplotlib": ("https://matplotlib.org/stable/", None),
    "sympy": ("https://docs.sympy.org/latest/", None),
    "scipy": ("https://docs.scipy.org/doc/scipy/", None),
}
intersphinx_disabled_domains = ["std"]

templates_path = ["_templates"]

pygments_style = "sphinx"

# -- Options for HTML output

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = "sphinx_book_theme"
html_logo = "gfx/logo.png"
html_favicon = ""
html_title = ""  # project
html_js_files = [
    "https://cdnjs.cloudflare.com/ajax/libs/require.js/2.3.4/require.min.js",
]
# sphinx-book-theme options
html_theme_options = {
    "repository_url": f"https://github.com/{github_user}/{github_repo_name}",
    "use_repository_button": True,
    "use_issues_button": True,
    "use_edit_page_button": True,
    "repository_branch": f"{github_version}",
    "path_to_docs": "docs",
    "use_download_button": True,
    "home_page_in_toc": True,
    "launch_buttons": {
        "binderhub_url": "https://mybinder.org",
        "thebe": True,
        "notebook_interface": "jupyterlab",
    },
}
html_sidebars = {
    "**": [
        "sidebar-logo.html",
        "searchbox.html",
        "sbt-sidebar-nav.html",
    ]
}
