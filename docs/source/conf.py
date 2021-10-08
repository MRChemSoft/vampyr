# Configuration file for the Sphinx documentation builder.

from pathlib import Path


def version_from_file(f):
    with f.open("r") as fp:
        return fp.read().strip()


# -- Project information

project = "VAMPyR"
copyright = "TO BE FILLED"
author = "TO BE FILLED"

version_file = Path(__file__).absolute().parents[2] / "VERSION"
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
]

autosummary_generate = True

intersphinx_mapping = {
    "python": ("https://docs.python.org/3/", None),
    "sphinx": ("https://www.sphinx-doc.org/en/master/", None),
}
intersphinx_disabled_domains = ["std"]

source_suffix = ".rst"

templates_path = ["_templates"]

pygments_style = "sphinx"

# -- Options for HTML output

html_theme = "sphinx_rtd_theme"

# -- Options for EPUB output
epub_show_urls = "footnote"
