function [K, geometry, info] = loadKernels(filename)
% LOADKERNELS  Read kernels written by saveKernels from an HDF5 file.
%
%   K = loadKernels(filename)
%   [K, geometry, info] = loadKernels(filename)
%
%   Returns the kernel struct exactly as it was saved (e.g. evl with KK,
%   KL, LK, LL and, if present, certificate and options; or devl with KO,
%   LO), with every array in its MATLAB shape and class.
%
% INPUTS
%   filename  path of an .h5 file written by saveKernels
%
% OUTPUTS
%   K         kernel struct (everything except /geometry)
%   geometry  struct from /geometry (fault, mantle, obs), or [] if the
%             file has none
%   info      struct of the root attributes (format, format_version,
%             index_order, created, description)
%
% SEE ALSO
%   saveKernels

filename = char(filename);
if ~isfile(filename)
    error('loadKernels:missing', 'Kernel file not found: %s', filename);
end
I = h5info(filename);
info = struct();
for a = I.Attributes'
    info.(a.Name) = a.Value;
end
if ~isfield(info, 'format') || ~strcmp(info.format, 'I2EC kernels')
    error('loadKernels:format', '%s was not written by saveKernels.', filename);
end

K = readGroup(filename, I, '/geometry');
geometry = [];
for g = I.Groups'
    if strcmp(g.Name, '/geometry')
        geometry = readGroup(filename, g, '');
    end
end
end

function S = readGroup(filename, G, skip)
S = struct();
for d = G.Datasets'
    path = [regexprep(G.Name, '/$', '') '/' d.Name];
    B = h5read(filename, path);
    sz = double(h5readatt(filename, path, 'matlab_size'))';
    if numel(sz) > 2 || all(sz > 1)
        A = permute(B, ndims(B):-1:1);
    else
        A = B;
    end
    A = reshape(A, sz);
    S.(d.Name) = A;
end
for g = G.Groups'
    if strcmp(g.Name, skip)
        continue
    end
    [~, name] = fileparts(g.Name);
    S.(name) = readGroup(filename, g, skip);
end
end
