# Example 3: mesh convergence

```matlab
run examples/example3_meshConvergence
```

The script solves the cycle of {doc}`example2` on three mantle meshes of the
same domain, made with `meshing/createShearZoneMesh` with largest edges of
40, 20 and 10 km (127, 495 and 1955 cells). The fault and boundary meshes are
the same, and the semi-analytical solver is exact in time, so any difference
between the three solutions comes from the mantle mesh. The script prints the
relative difference from the fine mesh at $t/T = 0.01$, $0.1$ and $0.99$ and
plots strain-rate snapshots and profiles on the three meshes.

The kernels of the three meshes are cached in `kernels/meshconvergence/`
(git-ignored). The first run computes them, about 15 minutes, mostly for the
fine mesh; later runs take about a minute.
