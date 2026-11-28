(kernels)=
# Kernels

The kernels are **half-space** solutions with a **boundary correction**:

- **Fault patches** use the plane-strain line-dislocation solutions of Davis &
  Martel.
- **Mantle cells** use the analytic stress of a uniformly strained triangle:
  the exact derivative of Barbot's (2018) closed-form displacement. It agrees
  with an independent quadrature to about $10^{-12}$ and satisfies Betti
  reciprocity to rounding.
- **The model box** is imposed by a boundary-element solve: zero displacement
  on its fixed sides and zero traction on its free sides. The background rate
  comes from the long-term solution, not from boundary motion.

:::{admonition} In the code
:class: in-code

{mat:func}`computeStressPlaneStrainTriangleShearZone` is the
triangle stress, generated with the Symbolic Math Toolbox by
`derivation/generateTriangleStress.m`; near a few lines through its edges,
where the closed form loses precision, that edge's contribution is integrated
numerically. {mat:func}`computeShzBoundaryCorrection` solves the boundary
problem once for all mantle sources. The fixed and free sides are set by
`boundary.Vpl` (0 and 1) in `examples/inputs/boundary2d.seg`.
:::

(mantle-kernels)=
## Mantle kernels: cell-averaged and certified

Physically, the stress caused by a viscous strain always opposes it, so
viscous relaxation can only decay. Sampling the stress at each cell's centre
breaks this between neighbouring cells and creates a few spurious growing
modes (43 on the example mesh), however exactly the stress is computed.

I2EC therefore **averages the stress over each cell**, a Galerkin
discretisation that keeps the right symmetry and sign, and **certifies** the
result: the mantle–mantle kernel is checked to have no growing mode and no
zero-energy (neutral) mode. A certified kernel is stable for any viscosity and
any power-law exponent, in both solvers.

:::{admonition} In the code
:class: in-code

this is the default of {mat:func}`computeAllStressKernelsBem`
(options `cellaverage` and `kernelcertify`, check by
{mat:func}`certifyShzKernel`, recorded in `evl.certificate`). Centre-sampled
kernels (`'cellaverage', 1, 'kernelcertify', false`) remain for comparison.
Their growing modes need a fix: `kernelmodify` ({mat:func}`modifyGrowingModes`)
for `ode45`, or {mat:func}`negateGrowingModes` on $A$, with every forcing term
built from the negated operator. Their results at the cell scale should not be
interpreted.
:::
