function m = loadExampleModel(withKernels)
%LOADEXAMPLEMODEL  The example geometry in examples/inputs/, as used by the MAIN scripts.
%
%   m = loadExampleModel()            geometry only
%   m = loadExampleModel(true)        also load examples/kernels/evl_orig.h5 and devl.h5
%
% Output struct m: em, rcv (49 patches), bnd (112), shz (521 triangles),
% and with kernels: evl (evl_orig), devl, obs (the 401 surface points that
% examples/exampleModel.m pairs with devl.h5).

if nargin < 1
    withKernels = false;
end

root = i2ecRoot();
m = struct();
m.em = geometry.LDhs(30e3, 0.25);
m.rcv = geometry.receiver(fullfile(root, 'examples', 'inputs', 'megathrust2d.seg'), m.em);
m.bnd = geometry.receiver(fullfile(root, 'examples', 'inputs', 'boundary2d.seg'), m.em);
m.shz = geometry.shearZoneReceiver(fullfile(root, 'examples', 'inputs', 'shearzone'), m.em);

if withKernels
    m.evl = loadKernels(fullfile(root, 'examples', 'kernels', 'evl_orig.h5'));
    m.devl = loadKernels(fullfile(root, 'examples', 'kernels', 'devl.h5'));
    m.obs = ([1; 0] * linspace(-100, 500, 401))' * 1e3;
end

end
