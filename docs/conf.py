# Sphinx configuration for the I2EC (2D) documentation.
#
# Build with docs/build.sh. Only library code is documented: the MATLAB
# source directory is docs/_build/matlab_src/, which holds symlinks to the
# code folders listed in MATLAB_CODE_DIRS (never data, scratch or example
# folders).
import os
import shutil
import sys
from pathlib import Path

DOCS = Path(__file__).resolve().parent
ROOT = DOCS.parent
sys.path.insert(0, str(DOCS / "_ext"))

# -- Project ------------------------------------------------------------------
project = "I2EC 2D"
author = "Rishav Mallick"
copyright = "2023–2026, Rishav Mallick"
release = ""
REPO_URL = "https://github.com/mallickrishg/I2EC_2d"

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.mathjax",
    "sphinxcontrib.matlab",
    "myst_parser",
    "matlab_help",          # docs/_ext/matlab_help.py: MATLAB help text -> reST
]

# -- MATLAB source (symlinks to the code folders only) -------------------------
MATLAB_CODE_DIRS = ["functions", "+geometry"]
MATLAB_SRC = DOCS / "_build" / "matlab_src"


def _link_code_dirs():
    MATLAB_SRC.mkdir(parents=True, exist_ok=True)
    for entry in MATLAB_SRC.iterdir():          # drop links to folders no longer listed
        if entry.name not in MATLAB_CODE_DIRS:
            entry.unlink() if entry.is_symlink() or entry.is_file() else shutil.rmtree(entry)
    for name in MATLAB_CODE_DIRS:
        link, target = MATLAB_SRC / name, ROOT / name
        if link.is_symlink() and Path(os.readlink(link)) != target:
            link.unlink()                         # points elsewhere: remake it
        elif link.exists() and not link.is_symlink():
            shutil.rmtree(link)
        if not link.is_symlink():
            link.symlink_to(target, target_is_directory=True)


_link_code_dirs()
matlab_src_dir = str(MATLAB_SRC)
primary_domain = "mat"
matlab_short_links = True
matlab_keep_package_prefix = False
matlab_auto_link = None
# parameter sizes such as "n x 1" are types, not objects to link to
nitpick_ignore_regex = [("mat:obj", r".*\bx\b.*")]

# -- MyST (Markdown pages) -----------------------------------------------------
myst_enable_extensions = ["dollarmath", "amsmath", "colon_fence"]
myst_heading_anchors = 3
suppress_warnings = ["myst.header"]

# -- General -------------------------------------------------------------------
source_suffix = {".rst": "restructuredtext", ".md": "markdown"}
templates_path = ["_templates"]
exclude_patterns = ["_build", "_ext", "figures", "dev", "README.md", "Thumbs.db", ".DS_Store"]

# -- HTML ----------------------------------------------------------------------
html_theme = "furo"
html_title = "I2EC 2D"
html_static_path = ["_static"]
html_css_files = ["custom.css"]
html_baseurl = "https://mallickrishg.github.io/I2EC_2d/"
html_theme_options = {
    "top_of_page_buttons": [],              # no "edit this page" / "view source"
    "footer_icons": [
        {
            "name": "GitHub",
            "url": REPO_URL,
            "html": (
                '<svg stroke="currentColor" fill="currentColor" stroke-width="0" '
                'viewBox="0 0 16 16"><path fill-rule="evenodd" d="M8 0C3.58 0 0 '
                "3.58 0 8c0 3.54 2.29 6.53 5.47 7.59.4.07.55-.17.55-.38 "
                "0-.19-.01-.82-.01-1.49-2.01.37-2.53-.49-2.69-.94-.09-.23-.48-.94-.82"
                "-1.13-.28-.15-.68-.52-.01-.53.63-.01 1.08.58 1.23.82.72 1.21 1.87.87 "
                "2.33.66.07-.52.28-.87.51-1.07-1.78-.2-3.64-.89-3.64-3.95 0-.87.31-1.59"
                ".82-2.15-.08-.2-.36-1.02.08-2.12 0 0 .67-.21 2.2.82.64-.18 1.32-.27 "
                "2-.27.68 0 1.36.09 2 .27 1.53-1.04 2.2-.82 2.2-.82.44 1.1.16 1.92.08 "
                "2.12.51.56.82 1.27.82 2.15 0 3.07-1.87 3.75-3.65 3.95.29.25.54.73.54 "
                "1.48 0 1.07-.01 1.93-.01 2.2 0 .21.15.46.55.38A8.013 8.013 0 0016 8c0"
                '-4.42-3.58-8-8-8z"></path></svg>'
            ),
            "class": "",
        },
    ],
}
