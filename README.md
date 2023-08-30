# README #

# Integral Imposed Earthquake Cycles for 2-d subduction zones

Authors: Rishav Mallick (Caltech Seismological Laboratory) & Sharadha Sathiakumar (Earth Observatory of Singapore)

MATLAB code to simulate post- and inter-seismic evolution of fault velocity and strain rate in a visco-elasto medium with a frictional fault. The long-term strain-rates need to be provided by geodynamic solutions or with known solutions to corner flow problems, while the short-term evolution is solved using an integral scheme. We use RK-4th order integration to solve the governing system of coupled ODEs (IVP).

This code uses code from the following repositories:

- a modification/version of mesh2d - https://github.com/dengwirda/mesh2d provided by https://people.sc.fsu.edu/~jburkardt/classes/dis_2014/mesh2d/mesh2d.html
- displacement and stress kernels for eigen sources in a linear elastic half-space (plane strain) from https://bitbucket.org/sbarbot/bssa-2018058/src/master/
- displacement and stress kernels for faults from https://github.com/Timmmdavis/CutAndDisplace
