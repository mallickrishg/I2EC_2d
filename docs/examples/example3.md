# Example 3: mesh convergence

The cycle of {doc}`example2` on three mantle meshes of the same domain, to see
how much the answer depends on the mesh. The semi-analytical solver is exact in
time and the fault and boundary are the same, so any difference comes from the
mantle mesh. Run it with

```matlab
run examples/example3_meshConvergence
```

The first run builds and caches the meshes and their kernels in
`examples/kernels/meshconvergence/` (git-ignored), which takes about a minute;
later runs load them. The steps below are the script, shortened.

## 1. Make the meshes

```matlab
res = struct('name', {'coarse', 'regular', 'fine'}, 'Hmax', {40, 20, 10});
% for each mesh k (steps 1 to 3 run in a loop over the three meshes):
name = sprintf('shearzone_h%02d', res(k).Hmax);
createShearZoneMesh('Hmax', res(k).Hmax, 'OutputFolder', cachedir, 'Name', name);
```

{mat:func}`createShearZoneMesh` (needs the PDE Toolbox)
triangulates the oceanic mantle and the mantle wedge with largest edges of 40,
20 and 10 km, and writes the `_vertices.dat` and `_triangulation.dat` files
described in {doc}`../inputs`.

## 2. Build or load the kernels of each mesh

```matlab
mk = exampleModel('ShearZone', fullfile(cachedir, name));
kernelfile = fullfile(cachedir, ['evl_' name '.h5']);
if isfile(kernelfile)
    evl = loadKernels(kernelfile);
else
    evl = computeAllStressKernelsBem(mk.rcv, mk.shz, mk.boundary);
    saveKernels(kernelfile, evl, 'rcv', mk.rcv, 'shz', mk.shz);
end
```

`exampleModel('ShearZone', ...)` loads the same fault and boundary with the new
mantle mesh. {mat:func}`computeAllStressKernelsBem` builds the default kernels
(cell-averaged, certified); {mat:func}`saveKernels` and
{mat:func}`loadKernels` store them as HDF5. This is the pattern to use for
your own geometry.

## 3. Solve the cycle on each mesh

```matlab
cyc = newtonianCycle(mk.rcv, mk.shz, evl, mk.Vpl, yr);
sol = cyc.solve(tfrac*cyc.Trecur);
res(k).edot = sqrt(sol.e22dot.^2 + sol.e23dot.^2);
```

The same rheology and solver as {doc}`example2`, at $t/T = 0.01$, $0.1$ and
$0.99$.

## 4. Compare the meshes

Each mesh is sampled on a common 2 km grid (the value of the triangle that
contains each point), and the coarse and regular solutions are compared with
the fine one:

```text
Difference from the fine mesh, 21605 points on a 2 km grid
                            t/T = 0.01    t/T = 0.10    t/T = 0.99
coarse   relative L2              28.6 %        54.2 %        57.7 %
         median |log10 ratio|        0.033         0.152         0.109
regular  relative L2              20.1 %        36.3 %        39.1 %
         median |log10 ratio|        0.016         0.068         0.046
```

The typical (median) difference roughly halves with each refinement. The L2
difference is dominated by the jumps between neighbouring cells where the
strain rate changes sharply, and decreases more slowly.

## 5. Plot the results

**Snapshots** on the three meshes (rows) at the three times (columns). The
large-scale pattern is the same on all three; the coarse mesh smears the
regions of fast and slow deformation.

```{figure} ../_static/examples/example3_meshConvergence_1.png
:alt: Mantle strain rate at three times on the coarse, regular and fine meshes
:width: 100%
```

**Profiles** at 60 km depth (top) and 300 km from the trench (bottom). The
regular and fine meshes agree closely; the coarse mesh departs most where the
strain rate changes sharply.

```{figure} ../_static/examples/example3_meshConvergence_2.png
:alt: Strain-rate profiles at 60 km depth and at 300 km from the trench on the three meshes
:width: 100%
```
