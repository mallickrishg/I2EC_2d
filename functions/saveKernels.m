function saveKernels(filename, K, varargin)
% SAVEKERNELS  Write kernels to an HDF5 file readable outside MATLAB.
%
%   saveKernels(filename, evl)
%   saveKernels(filename, devl)
%   saveKernels(filename, K, 'rcv', rcv, 'shz', shz, 'obs', obs, ...
%               'description', text)
%
%   Writes every numeric field of the struct K (e.g. KK, KL, LK, LL from
%   computeAllStressKernelsBem, or KO, LO from
%   computeAllDisplacementKernelsBem) as an HDF5 dataset of the same name.
%   Nested structs (evl.certificate, evl.options) become groups. Read the
%   file back with loadKernels.
%
%   Layout (format 'I2EC kernels', version 1):
%     /KK, /KL, ...          kernel blocks (float64, deflate-compressed)
%     /certificate/...       scalars from certifyMantleKernel, if present
%     /options/...           build options, if present
%     /geometry/fault/...    xc, W, dip of the fault patches  (if 'rcv')
%     /geometry/mantle/...   vertices, triangles, centroids, area  ('shz')
%     /geometry/obs          observation points  (if 'obs')
%   Every dataset has a 'description' and a 'units' attribute; the root
%   has 'format', 'format_version', 'index_order', 'created' and
%   'description'.
%
%   Index order: datasets are written so that a row-major (C-order)
%   reader such as h5py/NumPy sees the same shape and index order as
%   MATLAB: f['KL'][i, j, c] equals evl.KL(i+1, j+1, c+1). Triangle vertex
%   indices in /geometry/mantle/triangles are 1-based, as in MATLAB.
%
% INPUTS
%   filename  path of the .h5 file; an existing file is replaced
%   K         struct of kernel arrays
%
% OPTIONS
%   'rcv'          geometry.receiver of the fault, stored in /geometry
%   'shz'          geometry.shearZoneReceiver of the mantle cells
%   'obs'          [Nobs x 2] observation points [x, z] (m)
%   'description'  free text stored as a root attribute  ['']
%
% EXAMPLE
%   evl = computeAllStressKernelsBem(rcv, shz, boundary);
%   saveKernels('kernels/evl_orig.h5', evl, 'rcv', rcv, 'shz', shz);
%   % in Python:  import h5py; KK = h5py.File('kernels/evl_orig.h5')['KK'][:]
%
% SEE ALSO
%   loadKernels, computeAllStressKernelsBem, computeAllDisplacementKernelsBem

p = inputParser;
p.addParameter('rcv', []);
p.addParameter('shz', []);
p.addParameter('obs', [], @(x) isempty(x) || (isnumeric(x) && size(x, 2) == 2));
p.addParameter('description', '', @(s) ischar(s) || isstring(s));
p.parse(varargin{:});
opt = p.Results;

filename = char(filename);
if isfile(filename)
    delete(filename);
end

writeStruct(filename, '', K);

if ~isempty(opt.rcv)
    writeArray(filename, '/geometry/fault/xc', opt.rcv.xc, 'patch centres [x, z]', 'm');
    writeArray(filename, '/geometry/fault/W', opt.rcv.W, 'patch widths', 'm');
    writeArray(filename, '/geometry/fault/dip', opt.rcv.dip, 'patch dips', 'degrees');
end
if ~isempty(opt.shz)
    writeArray(filename, '/geometry/mantle/vertices', opt.shz.vert, 'mesh vertices [x, z]', 'm');
    writeArray(filename, '/geometry/mantle/triangles', int32(opt.shz.tri), ...
        'vertex indices of each triangle (1-based)', '');
    writeArray(filename, '/geometry/mantle/centroids', opt.shz.xc, 'cell centroids [x, z]', 'm');
    writeArray(filename, '/geometry/mantle/area', opt.shz.area, 'cell areas', 'm^2');
end
if ~isempty(opt.obs)
    writeArray(filename, '/geometry/obs', opt.obs, 'observation points [x, z]', 'm');
end

h5writeatt(filename, '/', 'format', 'I2EC kernels');
h5writeatt(filename, '/', 'format_version', int32(1));
h5writeatt(filename, '/', 'index_order', ['row-major: dataset[i, j, ...] = MATLAB array(i+1, j+1, ...); ' ...
    'h5py/NumPy see the MATLAB shape']);
h5writeatt(filename, '/', 'created', char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss')));
h5writeatt(filename, '/', 'description', char(opt.description));
end

function writeStruct(filename, prefix, S)
% numeric and logical fields -> datasets, struct fields -> groups
fields = fieldnames(S);
for k = 1:numel(fields)
    name = fields{k};
    value = S.(name);
    path = [prefix '/' name];
    if isstruct(value)
        writeStruct(filename, path, value);
    elseif isnumeric(value) || islogical(value)
        [descr, units] = describe(name);
        writeArray(filename, path, value, descr, units);
    else
        error('saveKernels:type', 'Field %s has class %s; only numeric, logical and struct fields are stored.', ...
            path(2:end), class(value));
    end
end
end

function writeArray(filename, path, A, descr, units)
if islogical(A)
    A = uint8(A);
end
sz = size(A);
if isvector(A) || isscalar(A)
    B = A(:);                              % 1-D dataset of length numel(A)
    dims = numel(A);
else
    B = permute(A, ndims(A):-1:1);         % C-order readers then see size(A)
    dims = size(B);
end
args = {'Datatype', class(B)};
if numel(B) > 1
    args = [args, {'ChunkSize', chunkSize(dims), 'Deflate', 4, 'Shuffle', true}]; %#ok<AGROW>
end
h5create(filename, path, dims, args{:});
h5write(filename, path, B);
h5writeatt(filename, path, 'description', descr);
h5writeatt(filename, path, 'units', units);
h5writeatt(filename, path, 'matlab_size', int64(sz));
end

function c = chunkSize(dims)
% chunks of at most 2^20 elements: halve the largest dimension
c = dims;
while prod(c) > 2^20
    [~, i] = max(c);
    c(i) = ceil(c(i)/2);
end
end

function [descr, units] = describe(name)
switch name
    case 'KK'
        descr = 'shear traction on fault patch i per unit slip on patch j';
        units = 'MPa/m';
    case 'KL'
        descr = ['deviatoric stress in mantle cell i per unit slip on fault patch j; ' ...
            'component c: 1 = (sxx - szz)/2, 2 = sxz'];
        units = 'MPa/m';
    case 'LK'
        descr = ['shear traction on fault patch i per unit deviatoric strain in mantle cell j; ' ...
            'component c: 1 = e22 (with e33 = -e22), 2 = e23'];
        units = 'MPa';
    case 'LL'
        descr = ['deviatoric stress in mantle cell i per unit deviatoric strain in cell j; ' ...
            'index 3 = stress component (1 = (sxx - szz)/2, 2 = sxz), index 4 = strain component (1 = e22, 2 = e23)'];
        units = 'MPa';
    case 'KO'
        descr = 'displacement at observation point i per unit slip on fault patch j; component c: 1 = ux, 2 = uz';
        units = 'm/m';
    case 'LO'
        descr = ['displacement at observation point i per unit deviatoric strain in mantle cell j; ' ...
            'index 3 = displacement component (1 = ux, 2 = uz), index 4 = strain component (1 = e22, 2 = e23)'];
        units = 'm';
    otherwise
        descr = name;
        units = '';
end
end
