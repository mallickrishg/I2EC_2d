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

## License

I2EC is released under the MIT License (see `LICENSE` in the repository).
