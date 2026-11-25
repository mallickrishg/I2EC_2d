Inverse problems
================

Gradients and sensitivities of the linear (Maxwell, Burgers) cycle with
respect to viscosity, for inversions of surface displacement (see
:doc:`../examples/example6`). These functions live in ``inverse/``. They
work in three layers: the forward model and its adjoint and tangent per
state row; the maps from rows to fault and mantle parameters; and the
station observations and misfit. The caller adds the last layer, the
parameterisation (e.g. low-rank KL modes), by the chain rule.

.. mat:module:: inverse

Forward, adjoint and tangent
----------------------------

.. mat:autofunction:: solveLinearSequence

.. mat:autofunction:: computeAdjointGradient

.. mat:autofunction:: computeTangentSequence

.. mat:autofunction:: computeExpDividedDifferences

Parameters
----------

.. mat:autofunction:: extractViscosityGradient

.. mat:autofunction:: assembleRowPerturbation

Observations
------------

.. mat:autofunction:: assembleStationKernels

.. mat:autofunction:: computeDisplacementMisfit
