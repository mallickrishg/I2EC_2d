function m = exampleModel(varargin)
% EXAMPLEMODEL  Geometry, meshes and kernels shared by the I2EC examples.
%
%   m = exampleModel()
%   m = exampleModel('Name', Value, ...)
%
%   Adds the I2EC folders to the path and loads the example subduction
%   zone: the megathrust (inputs/megathrust2d.seg), the model boundary
%   (inputs/boundary2d.seg), the kinematic hinge fault (inputs/hinge2d.seg)
%   and a mantle mesh (default inputs/shearzone). With 'StressKernels' or
%   'DisplacementKernels' true it also returns the default kernels: loaded
%   from the cache in kernels/ when they match the meshes, computed (and
%   not saved) otherwise.
%
% OPTIONS (name-value)
%   'ShearZone'            mesh prefix, relative to the repository root
%                          ['inputs/shearzone']
%   'StressKernels'        return m.evl (computeAllStressKernelsBem,
%                          default options)  [false]
%   'DisplacementKernels'  return m.devl and the hinge kernels m.Ghx,
%                          m.Ghz at the surface points m.obs  [false]
%
% OUTPUTS
%   m   struct with fields
%       .root        repository root
%       .mu, .nu     shear modulus (MPa) and Poisson's ratio
%       .Vpl         plate rate (m/s)
%       .yr          seconds per year (3.15e7, as in the rest of I2EC)
%       .earthModel  geometry.LDhs object
%       .rcv         geometry.receiver: megathrust (49 patches)
%       .boundary    geometry.receiver: model boundary
%       .hinge       geometry.receiver: kinematic hinge fault
%       .shz         geometry.shearZoneReceiver: mantle mesh
%       .obs         [401 x 2] surface points, x from -100 to 500 km (m)
%       .evl, .devl, .Ghx, .Ghz  kernels, when requested
%
% NOTES
%   - The cached kernels are kernels/evl_orig.mat (stress, default options)
%     and kernels/devl.mat (displacement at m.obs), both for
%     inputs/shearzone. They are used only when their sizes match the
%     meshes; a different mantle mesh gets its kernels computed.
%   - rcv, shz and the other receivers are handle objects: examples set
%     their rheology on the returned objects.
%
% SEE ALSO
%   newtonianCycle, computeAllStressKernelsBem,
%   computeAllDisplacementKernelsBem

p = inputParser;
p.addParameter('ShearZone', fullfile('inputs', 'shearzone'));
p.addParameter('StressKernels', false);
p.addParameter('DisplacementKernels', false);
p.parse(varargin{:});
opt = p.Results;

m.root = fileparts(fileparts(mfilename('fullpath')));
addpath(m.root, fullfile(m.root, 'functions'), fullfile(m.root, 'meshing'), ...
    fullfile(m.root, 'examples'));

m.mu = 30e3;      % MPa
m.nu = 0.25;
m.Vpl = 1e-9;     % m/s
m.yr = 3.15e7;    % s

m.earthModel = geometry.LDhs(m.mu, m.nu);
m.rcv = geometry.receiver(fullfile(m.root, 'inputs', 'megathrust2d.seg'), m.earthModel);
% Vpl column of the boundary file = boundary-condition type (0: zero
% displacement, 1: zero traction); the kernels solve the homogeneous problem
m.boundary = geometry.receiver(fullfile(m.root, 'inputs', 'boundary2d.seg'), m.earthModel);
m.hinge = geometry.receiver(fullfile(m.root, 'inputs', 'hinge2d.seg'), m.earthModel);
shzPrefix = opt.ShearZone;
if ~isfile([shzPrefix '_vertices.dat'])
    shzPrefix = fullfile(m.root, shzPrefix);
end
m.shz = geometry.shearZoneReceiver(shzPrefix, m.earthModel);
m.obs = [linspace(-100, 500, 401)' * 1e3, zeros(401, 1)];

if opt.StressKernels
    cached = fullfile(m.root, 'kernels', 'evl_orig.mat');
    if isfile(cached)
        S = load(cached, 'evl_orig');
    end
    if isfile(cached) && size(S.evl_orig.KK, 1) == m.rcv.N && size(S.evl_orig.LL, 1) == m.shz.N
        m.evl = S.evl_orig;
    else
        m.evl = computeAllStressKernelsBem(m.rcv, m.shz, m.boundary);
    end
end

if opt.DisplacementKernels
    cached = fullfile(m.root, 'kernels', 'devl.mat');
    if isfile(cached)
        S = load(cached, 'devl');
    end
    if isfile(cached) && isequal(size(S.devl.KO), [size(m.obs, 1) m.rcv.N 2]) ...
            && size(S.devl.LO, 2) == m.shz.N
        m.devl = S.devl;
    else
        m.devl = computeAllDisplacementKernelsBem(m.obs, m.rcv, m.shz, m.boundary, 1);
    end
    [m.Ghx, m.Ghz] = computeFaultDisplacementKernelsBem(m.hinge, m.obs, m.boundary, 1);
end
end
