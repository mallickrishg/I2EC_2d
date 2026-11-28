# Loading and governing equations

The earthquake cycle is a departure from steady plate motion. This page
defines that steady state and the equations for the departure.

## Long-term loading

Deformation is split into a **long-term** part and an **earthquake-cycle**
part. Over geological time the plates move steadily: the fault slips at the
plate rate $v^{\infty}$ and the mantle flows as corner flow (Newtonian, or
power-law with $n = 3$). The cycle is a time-dependent departure from this
state, driven by locking and imposed earthquakes.

:::{admonition} In the code
:class: in-code

{mat:func}`computeLongTermStrainRates` gives the corner-flow
strain rates in each cell, {mat:func}`computeLongTermVelocity` the velocity
field ({doc}`../examples/index` plots it).
:::

## Governing equations

The fault is a set of patches on a boundary $\Gamma$ and the mantle a set of
triangular cells $\Omega$. Slip and viscous strain that depart from the
long-term state change the stress everywhere:

$$
\dot{\sigma}_i(\mathbf{x}) =
\int_{\Gamma} K_i(\mathbf{x};\xi)\,\bigl(v(\xi) - v^{\infty}\bigr)\,d\xi
+ \iint_{\Omega} K_{ji}(\mathbf{x};\alpha)\,
  \bigl(\dot{\epsilon}_{\mathrm{v}_j}(\alpha) -
        \dot{\epsilon}^{\infty}_{\mathrm{v}_j}(\alpha)\bigr)\,d\alpha .
$$

The kernels $K$ are elastic Green's functions. Mantle flow is incompressible,
so each cell has two unknowns, the deviatoric strain rates
$\dot\epsilon_{22}$ and $\dot\epsilon_{23}$. Discretised, the kernels are four
matrices ($N$ fault patches, $M$ mantle cells):

| Kernel | Maps | Size |
| --- | --- | --- |
| `KK` | fault slip → fault shear traction | $N \times N$ |
| `KL` | fault slip → mantle stress | $M \times N \times 2$ |
| `LK` | mantle strain → fault shear traction | $N \times M \times 2$ |
| `LL` | mantle strain → mantle stress | $M \times M \times 2 \times 2$ |

:::{admonition} In the code
:class: in-code

{mat:func}`computeAllStressKernelsBem` returns all four, in MPa
per metre of slip or per unit strain. Mantle stresses are the deviatoric pair
$(\sigma_{22}-\sigma_{33})/2$ and $\sigma_{23}$.
:::
