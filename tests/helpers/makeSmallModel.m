function m = makeSmallModel(folder)
%MAKESMALLMODEL  Small synthetic subduction geometry for fast tests.
%
%   m = makeSmallModel(folder)
%
% Writes a 10-patch megathrust, a coarse 20-patch box boundary and a
% 24-triangle mantle shear zone into FOLDER (use a temporary folder) and
% loads them with the +geometry classes.  The geometry follows the example
% in inputs/ (same dip, same boundary box), only coarser, so every kernel
% routine runs in about a second.
%
% Output struct m:
%   em          geometry.LDhs earth model (G = 30e3 MPa, nu = 0.25)
%   rcv         geometry.receiver   megathrust (N = 10)
%   bnd         geometry.receiver   boundary   (N = 20, Dirichlet u = 0)
%   shz         geometry.shearZoneReceiver mantle (N = 24)
%   obs         [7 x 2] surface observation points
%   folder      where the input files were written

G = 30e3;
nu = 0.25;

faultFile = fullfile(folder, 'fault.seg');
fid = fopen(faultFile, 'w');
fprintf(fid, '# n Vpl x z Width Dip W0 qW\n');
fprintf(fid, '1 1.0 0.0 1.0 103077.6 14.0362 10307.76 1 0 0\n');
fclose(fid);

bndFile = fullfile(folder, 'bnd.seg');
fid = fopen(bndFile, 'w');
fprintf(fid, '# n Vpl x z Width Dip W0 qW bcx bcz\n');
fprintf(fid, '1 0.0 -500e3 0     500e3   90 100e3 1 0 0\n');
fprintf(fid, '2 0.0 -500e3 500e3 1250e3   0 125e3 1 0 0\n');
fprintf(fid, '3 0.0  750e3 500e3 500e3  -90 100e3 1 0 0\n');
fclose(fid);

[xg, zg] = meshgrid(linspace(60e3, 160e3, 5), linspace(-60e3, -120e3, 4));
vert = [xg(:) zg(:)];
tri = delaunay(vert(:,1), vert(:,2));
writematrix(vert, fullfile(folder, 'mantle_vertices.dat'));
writematrix(tri, fullfile(folder, 'mantle_triangulation.dat'));

m = struct();
m.folder = folder;
m.em = geometry.LDhs(G, nu);
m.rcv = geometry.receiver(faultFile, m.em);
m.bnd = geometry.receiver(bndFile, m.em);
m.shz = geometry.shearZoneReceiver(fullfile(folder, 'mantle'), m.em);
m.obs = [linspace(-50e3, 200e3, 7)' zeros(7,1)];

end
