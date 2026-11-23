function [vert, tri] = createShearZoneMesh(varargin)
% CREATESHEARZONEMESH  Mesh the oceanic mantle and mantle wedge into triangles.
%
%   createShearZoneMesh()
%   createShearZoneMesh('Name', Value, ...)
%   [vert, tri] = createShearZoneMesh(...)
%
%   Meshes two polygons with the MATLAB PDE Toolbox (decsg,
%   geometryFromEdges, generateMesh, linear triangles): the oceanic mantle
%   below the incoming plate and the slab, and the mantle wedge between the
%   upper plate and the slab. The two meshes are joined (wedge vertex
%   indices offset by the number of oceanic vertices) and written as
%       <OutputFolder>/<Name>_vertices.dat        x, z (m), one row per vertex
%       <OutputFolder>/<Name>_triangulation.dat   3 vertex indices per row
%   which geometry.shearZoneReceiver reads.
%
% OPTIONS (name-value)
%   'Hmax'          maximum edge length (km); scalar, or [oceanic, wedge].
%                   Default 20
%   'Oceanic'       [N x 2] oceanic-mantle polygon [x, z] (km)
%   'Wedge'         [N x 2] mantle-wedge polygon [x, z] (km)
%   'OutputFolder'  folder for the .dat files. Default: examples/inputs/ of the
%                   repository. Pass '' to write nothing
%   'Name'          file name prefix. Default 'shearzone'
%   'Plot'          true to plot the mesh. Default false
%
% OUTPUTS
%   vert            [P x 2] vertex coordinates [x, z] (m)
%   tri             [M x 3] vertex indices of each triangle
%
% NOTES
%   - The default polygons are the outline of the example mesh in examples/inputs/
%     (slab dip atan(0.25) = 14.04 deg, oceanic plate 20 km and upper
%     plate 35 km thick, domain to x = 450 km and z = -200 km). The wedge
%     starts at (200, -50) km, the down-dip end of the example megathrust,
%     so the fault ends on the corner vertex of the wedge.
%   - The example mesh in examples/inputs/ was made with the PDE Toolbox, but the
%     exact settings were not recorded: 'Hmax' = 20 gives 497 triangles,
%     close to its 521 but not the same mesh. Rebuild the cached kernels
%     in examples/kernels/ after replacing it.
%   - Requires the PDE Toolbox. It is needed only to make meshes, not to
%     run models.
%
% EXAMPLE
%   % finer mesh written next to the example one
%   createShearZoneMesh('Hmax', [15 10], 'Name', 'shearzone_fine', 'Plot', true);
%   shz = geometry.shearZoneReceiver('examples/inputs/shearzone_fine', geometry.LDhs(30e3, 0.25));
%
% SEE ALSO
%   geometry.shearZoneReceiver, generateMesh

p = inputParser;
p.addParameter('Hmax', 20, @(h) isnumeric(h) && any(numel(h) == [1 2]) && all(h > 0));
p.addParameter('Oceanic', [-100 -20; -100 -200; 450 -200; 450 -132.5; 0 -20], @isPolygon);
p.addParameter('Wedge', [200 -50; 450 -112.5; 450 -35; 200 -35], @isPolygon);
p.addParameter('OutputFolder', fullfile(fileparts(fileparts(mfilename('fullpath'))), 'examples', 'inputs'), ...
    @(s) ischar(s) || isstring(s));
p.addParameter('Name', 'shearzone', @(s) ischar(s) || isstring(s));
p.addParameter('Plot', false, @islogical);
p.parse(varargin{:});
opt = p.Results;

Hmax = opt.Hmax .* [1 1];

[vo, to] = meshPolygon(opt.Oceanic, Hmax(1));
[vw, tw] = meshPolygon(opt.Wedge, Hmax(2));

vert = [vo; vw] * 1e3;                       % km -> m
tri = [to; tw + size(vo, 1)];

if strlength(opt.OutputFolder) > 0
    base = fullfile(opt.OutputFolder, opt.Name);
    writematrix(vert, [char(base) '_vertices.dat']);
    writematrix(tri, [char(base) '_triangulation.dat']);
end

if opt.Plot
    figure
    triplot(tri, vert(:,1)/1e3, vert(:,2)/1e3, 'k-');
    axis equal tight
    xlabel('x (km)'); ylabel('z (km)');
    title(sprintf('%d triangles, %d vertices', size(tri, 1), size(vert, 1)));
end

end

function [v, t] = meshPolygon(P, Hmax)
model = createpde;
geometryFromEdges(model, decsg([2; size(P, 1); P(:, 1); P(:, 2)]));
mesh = generateMesh(model, 'Hmax', Hmax, 'GeometricOrder', 'linear');
v = mesh.Nodes';
t = double(mesh.Elements');
end

function tf = isPolygon(P)
tf = isnumeric(P) && size(P, 2) == 2 && size(P, 1) >= 3;
end
