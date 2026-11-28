function H = assembleStationKernels(devl, sys, varargin)
% ASSEMBLESTATIONKERNELS  Surface-displacement operator acting on the linear state.
%
%   H = assembleStationKernels(devl, sys)
%   H = assembleStationKernels(devl, sys, 'Components', 1, 'Stations', idx)
%
%   Builds H such that H*U is the surface displacement at the stations
%   due to the integrated state U (slip on creeping patches and strain,
%   e.g. fwd.U of solveLinearSequence). Kelvin rows add to the total
%   strain (H = G*sys.P), so it works for Maxwell and Burgers states.
%   Rows are ordered component-major: all stations of the first
%   component, then all of the next.
%
% INPUTS
%   devl   displacement kernels of computeAllDisplacementKernelsBem
%          (.KO [Nobs x N x 2], .LO [Nobs x M x 2 x 2])
%   sys    struct from assembleLinearSystem
%   'Components'  displacement components, 1 = x, 2 = z  [[1 2]]
%   'Stations'    station indices into the rows of devl  [all]
%
% OUTPUTS
%   H      [numel(Stations)*numel(Components) x n] operator
%
% NOTES
%   - Only the creeping patches are in the state; displacement from slip
%     on locked patches (coseismic slip, slip deficit) does not depend on
%     viscosity and belongs in the fixed part of computeDisplacementMisfit.
%
% SEE ALSO
%   computeAllDisplacementKernelsBem, computeDisplacementMisfit,
%   solveLinearSequence
%
% AUTHORS
%   Rishav Mallick (EOS), 2026

ip = inputParser;
ip.addParameter('Components', [1 2]);
ip.addParameter('Stations', 1:size(devl.KO, 1));
ip.parse(varargin{:});
comp = ip.Results.Components; sta = ip.Results.Stations;

G = [];
for c = comp(:)'
    G = [G; devl.KO(sta, ~sys.locked, c), devl.LO(sta, :, c, 1), devl.LO(sta, :, c, 2)]; %#ok<AGROW>
end
H = G * sys.P;
end
