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
(git-ignored). The first run computes them in about a minute (most of it for
the fine mesh); later runs load them and take about a minute.

The script prints, for points on a 2 km grid inside the mantle, the
difference of the coarse and regular solutions from the fine one:

| Mesh | | $t/T = 0.01$ | $t/T = 0.1$ | $t/T = 0.99$ |
| --- | --- | --- | --- | --- |
| coarse (127 cells) | relative L2 | 28% | 54% | 58% |
| | median $\lvert\log_{10}\text{ratio}\rvert$ | 0.033 | 0.151 | 0.107 |
| regular (495 cells) | relative L2 | 21% | 39% | 39% |
| | median $\lvert\log_{10}\text{ratio}\rvert$ | 0.016 | 0.068 | 0.045 |

The typical (median) difference roughly halves with each refinement: on the
regular mesh it is 4% just after the earthquake and about 10–17% later in the
cycle. The L2 difference is dominated by the jumps between neighbouring cells
where the strain rate changes sharply, and decreases more slowly.

```{figure} ../_static/examples/example3_meshConvergence_1.png
:alt: Mantle strain rate at three times on the coarse, regular and fine meshes
:width: 100%

Mantle strain rate at $t/T = 0.01$, $0.1$ and $0.99$ (columns) on the coarse,
regular and fine meshes (rows). The large-scale pattern is the same on all
three; the coarse mesh smears the regions of fast and slow deformation.
```

```{figure} ../_static/examples/example3_meshConvergence_2.png
:alt: Strain-rate profiles at 60 km depth and at x = 300 km on the three meshes
:width: 100%

Strain-rate profiles along $z = -60$ km (top) and $x = 300$ km (bottom) at the
same three times. The regular and fine meshes agree closely; the coarse mesh
departs most where the strain rate changes sharply.
```
