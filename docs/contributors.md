# Contributors

- **Rishav Mallick** built and maintains I2EC.
- **Sharadha Sathiakumar** co-developed the shear-zone classes, the long-term
  strain-rate and velocity calculations and the surface-displacement kernels,
  and the negated-mode exact solver and its tests.

## How to cite

If you use this code in your research, please cite:

> Mallick, R., & Sathiakumar, S. (2026). A fast numerical routine to model
> viscoelastic earthquake cycles in subduction zones. *Earth, Planets and
> Space*. <https://doi.org/10.1186/s40623-026-02429-4>

## Third-party components

This repository integrates and modifies components from:

- [MATLAB PDE Toolbox](https://www.mathworks.com/products/pde.html): 2D geometry
  discretisation and mesh generation.
- Barbot (2018): displacement and stress kernels for eigen sources in a linear
  elastic half-space (plane strain), <https://doi.org/10.1785/0120180058>.
- [CutAndDisplace](https://github.com/Timmmdavis/CutAndDisplace): kernels for
  fault displacement and stress.
- [bluewhitered](https://www.mathworks.com/matlabcentral/fileexchange/4058-bluewhitered)
  by Nathan Childress (MATLAB File Exchange, 2003): the blue-white-red colour
  map.
- [cline](https://www.mathworks.com/matlabcentral/fileexchange/14677-cline) by
  Sebastian Hölz (MATLAB File Exchange, 2007): lines coloured by value, used by
  {mat:func}`plotPatches2d`.
- `drchrnd`: the standard construction of a Dirichlet sample from normalised
  gamma variates (needs the Statistics and Machine Learning Toolbox).

(code-conventions)=
## Code conventions

For contributors adding or renaming functions:

1. **Functions** are lowerCamelCase and start with a verb: `compute…`,
   `assemble…`, `solve…`, `run…`, `integrate…`, `extract…`, `plot…`,
   `load…`/`save…`, `create…`, `certify…`, `is…`. No snake_case and no
   `get…`.
2. **Acronyms** are one capitalised word (`Bem`, `Shz`, `Ivp`, `Ode`), except
   long-established ones (`LD…`, `KK`/`KL`/`LK`/`LL`).
3. **Half-space kernels** live in `+geometry/` with no suffix; their
   boundary-corrected versions live in `functions/` with the suffix `Bem`, for
   example `geometry.computeFaultTractionKernels` and
   {mat:func}`computeFaultTractionKernelsBem`.
4. **Domain words:** `Fault` for megathrust patches, `Shz` for mantle cells
   (as the `shz` variable), `Boundary` for the model edge.
5. **Classes** keep their lowercase names (`geometry.receiver`,
   `shearZoneReceiver`, `earthModel`, …): every script uses them, and a
   case-only rename cannot keep the old name working on macOS or Windows.
6. **Third-party files** keep their names (`LDdispHS`, `LDstressHS`, `cline`,
   `bluewhitered`, `drchrnd`), as does `filenameGenerator`.
7. **Scripts** are named `exampleN_whatItShows` or `plot…`.

Functions renamed to follow these rules in November 2026 are listed in
{doc}`migration`.

## License

I2EC is released under the MIT License (see `LICENSE` in the repository).
