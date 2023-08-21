# README #

# Integral Imposed Earthquake Cycles for 2-d subduction zones

MATLAB code to simulate post- and inter-seismic evolution of fault velocity and strain rate in a visco-elasto medium with a frictional fault. The long-term strain-rates need to be provided by geodynamic solutions or with known solutions to corner flow problems, while the short-term evolution is solved using an integral scheme. We use RK-4th order integration to solve the governing system of coupled ODEs (IVP).

This code relies on 

a modification/version of mesh2d - https://github.com/dengwirda/mesh2d provided by https://people.sc.fsu.edu/~jburkardt/classes/dis_2014/mesh2d/mesh2d.html. This is incorporated into the package.
