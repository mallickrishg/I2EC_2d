# Prepare the model geometry

Example inputs are in `inputs/`. Lengths are in metres.

## Faults and boundaries (`.seg` files)

Faults and model boundaries are described as straight segments, each split into
patches. Read them with {mat:class}`geometry.receiver <receiver>`:

```matlab
earthModel = geometry.LDhs(mu, nu);
rcv      = geometry.receiver('inputs/megathrust2d.seg', earthModel);
boundary = geometry.receiver('inputs/boundary2d.seg',   earthModel);
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
example, `inputs/boundary2d.seg`:

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
shz = geometry.shearZoneReceiver('inputs/shearzone', earthModel);
```

New meshes can be made with `meshing/CREATE_shearzone_mesh.m` or
`meshing/CREATE_geometryFromEdgesMesh.m`.

## Pre-computed kernels

`kernels/evl_orig.mat` and `kernels/devl.mat` hold the stress and displacement
kernels for the example geometry. They are valid only for these meshes:
after any change in geometry, recompute them with
{mat:func}`computeAllStressKernelsBem` and
{mat:func}`computeAllDisplacementKernelsBem`.
