# Prepare the model geometry

Example inputs are in `examples/inputs/`. Lengths are in metres.

## Faults and boundaries (`.seg` files)

Faults and model boundaries are described as straight segments, each split into
patches. Read them with {mat:class}`geometry.receiver <receiver>`:

```matlab
earthModel = geometry.LDhs(mu, nu);
rcv      = geometry.receiver('examples/inputs/megathrust2d.seg', earthModel);
boundary = geometry.receiver('examples/inputs/boundary2d.seg',   earthModel);
```

Lines starting with `#` are comments. Each data line has ten columns:

| Column | Meaning |
| --- | --- |
| `n` | Segment index. |
| `Vpl` | Long-term slip rate of the segment (scaled by the plate rate in the driver scripts). |
| `x`, `z` | Upper end of the segment (m); `z` is depth, positive downward. |
| `Width` | Down-dip length of the segment (m). |
| `Dip` | Dip angle (degrees). |
| `W0` | Approximate width of the first patch (m). |
| `qW` | Geometric growth factor of patch width with distance (1 = uniform). |
| `bcx`, `bcz` | Boundary-condition velocities in *x* and *z* (for the model boundary; scaled by the plate rate in the driver scripts). |

{mat:func}`geometry.flt2flt <flt2flt>` splits each segment into patches. For
example, `examples/inputs/boundary2d.seg`:

```text
# n  Vpl     x       z     Width    Dip    W0   qW    bcx      bcz
1    0.0  -500.0e3  0.0e3   500e3     90   20e3   1    0         0
2    0.0  -500.0e3  500e3  1250e3      0   20e3   1    0         0
3    0.0   750.0e3  500e3   500e3    -90   20e3   1    0         0
```

To build a planar megathrust in MATLAB instead of a file, see
{mat:func}`create_megathrust`.

## Mantle mesh (`.dat` files)

The mantle is a triangular mesh stored in two comma-separated files that share
a prefix:

- `<name>_vertices.dat`: one `x,z` pair per vertex (m; `z` negative below the
  surface);
- `<name>_triangulation.dat`: three 1-based vertex indices per triangle.

Load them by passing the prefix to
{mat:class}`geometry.shearZoneReceiver <shearZoneReceiver>`:

```matlab
shz = geometry.shearZoneReceiver('examples/inputs/shearzone', earthModel);
```

New meshes can be made with {mat:func}`createShearZoneMesh`, which meshes the
oceanic mantle and the mantle wedge with the PDE Toolbox and writes both files:

```matlab
addpath functions
createShearZoneMesh('Hmax', [15 10], 'Name', 'shearzone_fine');   % examples/inputs/shearzone_fine_*.dat
```

Its default polygons are the outline of the example mesh. The example mesh was
made with the PDE Toolbox, but its exact settings were not recorded, so
`'Hmax', 20` gives a similar mesh (495 triangles instead of 521), not the same
one. Rebuild the cached kernels below after changing the mesh.

## Pre-computed kernels

`examples/kernels/evl_orig.h5` and `examples/kernels/devl.h5` hold the stress and displacement
kernels for the example geometry. `evl_orig.h5` holds the default kernels of
{mat:func}`computeAllStressKernelsBem` (3-point cell average, certified; see its
`certificate` group). They are valid only for these meshes:
after any change in geometry, recompute them with
{mat:func}`computeAllStressKernelsBem` and
{mat:func}`computeAllDisplacementKernelsBem`.

(kernel-options)=
### Kernel options

{mat:func}`computeAllStressKernelsBem` computes the mantle stress from the
analytic triangle solution and, by default, averages it over each cell and
certifies it (see {ref}`mantle-kernels`):

| Call | Mantle stress | Growing modes | Use |
| --- | --- | --- | --- |
| `computeAllStressKernelsBem(rcv, shz, bnd)` | averaged over each cell (3 points), certified | none | default, both solvers |
| `..., 'cellaverage', 7` | averaged (7 points), certified | none | finer quadrature, about twice the cost |
| `..., 'cellaverage', 1, 'kernelcertify', false` | at the cell centre | a few | comparison only |

Centre-sampled kernels need a fix for their growing modes: `'kernelmodify', 1`
or `2` for `ode45`, or `[~, ~, A] = negateGrowingModes(A)` for the
semi-analytical solvers.

### Kernel files (HDF5)

Kernels are stored as HDF5 files, so they can be read outside MATLAB.
{mat:func}`saveKernels` writes them and {mat:func}`loadKernels` reads them back
into the same MATLAB struct:

```matlab
evl = computeAllStressKernelsBem(rcv, shz, boundary);
saveKernels('examples/kernels/evl_orig.h5', evl, 'rcv', rcv, 'shz', shz);
evl = loadKernels('examples/kernels/evl_orig.h5');
```

Layout (format `I2EC kernels`, version 1):

| Path | Contents |
| --- | --- |
| `/KK`, `/KL`, `/LK`, `/LL` | stress kernel blocks (from {mat:func}`computeAllStressKernelsBem`) |
| `/KO`, `/LO` | displacement kernel blocks (from {mat:func}`computeAllDisplacementKernelsBem`) |
| `/certificate/...` | the check by {mat:func}`certifyMantleKernel`, if the kernel was certified |
| `/options/...` | `cellaverage`, `kernelcertify`, `kernelmodify` used to build it |
| `/geometry/fault/xc`, `W`, `dip` | fault patch centres (m), widths (m), dips (degrees) |
| `/geometry/mantle/vertices`, `triangles`, `centroids`, `area` | mantle mesh; triangle indices are 1-based |
| `/geometry/obs` | observation points of the displacement kernels (m) |

Every dataset has `description` and `units` attributes. Arrays are written so
that a row-major reader sees the MATLAB shape and index order:
`f['KL'][i, j, c]` is MATLAB's `evl.KL(i+1, j+1, c+1)`. Scalars (the
certificate values and options) are one-element arrays. In Python:

```python
import h5py
with h5py.File("examples/kernels/evl_orig.h5", "r") as f:
    KK = f["KK"][:]                      # (N, N)
    LL = f["LL"][:]                      # (M, M, 2, 2)
    print(f["LL"].attrs["description"], f["LL"].attrs["units"])
    area = f["geometry/mantle/area"][:]
```
