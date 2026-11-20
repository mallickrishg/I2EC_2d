function m = loadExampleModel(withKernels)
%LOADEXAMPLEMODEL  The example geometry in inputs/, as used by the MAIN scripts.
%
%   m = loadExampleModel()            geometry only
%   m = loadExampleModel(true)        also load kernels/evl_orig.mat and devl.mat
%
% Output struct m: em, rcv (49 patches), bnd (112), shz (521 triangles),
% and with kernels: evl (evl_orig), devl, obs (the 401 surface points that
% examples/exampleModel.m pairs with devl.mat).

if nargin < 1
    withKernels = false;
end

root = i2ecRoot();
m = struct();
m.em = geometry.LDhs(30e3, 0.25);
m.rcv = geometry.receiver(fullfile(root, 'inputs', 'megathrust2d.seg'), m.em);
m.bnd = geometry.receiver(fullfile(root, 'inputs', 'boundary2d.seg'), m.em);
m.shz = geometry.shearZoneReceiver(fullfile(root, 'inputs', 'shearzone'), m.em);

if withKernels
    S = load(fullfile(root, 'kernels', 'evl_orig.mat'), 'evl_orig');
    m.evl = S.evl_orig;
    S = load(fullfile(root, 'kernels', 'devl.mat'), 'devl');
    m.devl = S.devl;
    m.obs = ([1; 0] * linspace(-100, 500, 401))' * 1e3;
end

end
