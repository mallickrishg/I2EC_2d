# Build the documentation

The documentation uses [Sphinx](https://www.sphinx-doc.org) with the
[Furo](https://pradyunsg.me/furo/) theme,
[sphinxcontrib-matlabdomain](https://github.com/sphinx-contrib/matlabdomain) for
the MATLAB API and [MyST](https://myst-parser.readthedocs.io) for the Markdown
pages.

## Set up

Create a virtual environment outside the repository and install the
requirements:

```bash
python3 -m venv ~/.venvs/matlab-docs
```

```bash
~/.venvs/matlab-docs/bin/pip install -r docs/requirements.txt
```

## Build

```bash
docs/build.sh
```

The HTML goes to `docs/_build/html`. Add `--open` to open it in a browser, or
`--serve` to serve it at <http://127.0.0.1:8000>. The script uses
`~/.venvs/matlab-docs` if it exists, `$SPHINXBUILD` if set, and otherwise
`sphinx-build` on the `PATH`.

## How the API pages are made

- Only library code is documented. `docs/conf.py` fills
  `docs/_build/matlab_src/` with links to `functions/` and `+geometry/` and
  points the MATLAB domain there.
- The help text at the top of each `.m` file is converted to reStructuredText
  by `docs/_ext/matlab_help.py`. It recognises the usual MATLAB layout: a
  summary line, usage lines (dropped), `INPUTS`/`OUTPUTS` lists (shown as
  parameters and returns, with `[n x 1]` sizes as types), and `NOTES`,
  `EXAMPLE`, `SEE ALSO` and `AUTHOR` sections. The `.m` files are never
  changed.
- To check the conversion of every function, run
  `~/.venvs/matlab-docs/bin/python docs/_ext/check_help.py` (add
  `--show NAME` or `--show all` to print the generated reStructuredText).

## Figures

Figures are made by MATLAB scripts in `docs/figures/make_*.m`, which write PNG
files to `docs/_static/`. Rerun a script after changing the inputs it draws:

```matlab
run docs/figures/make_overview.m
```

## Publishing

The GitHub Actions workflow `.github/workflows/docs.yml` builds the pages and
publishes them to <https://mallickrishg.github.io/I2EC_2d/> on every push to
`main`.
