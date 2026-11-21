function rcv = createMegathrust(earthModel,x0,z0,dip,Fault_width,Nfault)
% CREATEMEGATHRUST  Create a planar megathrust and load it as patches.
%
%   rcv = createMegathrust(earthModel, x0, z0, dip, Fault_width, Nfault)
%
%   Writes a one-segment fault file to a temporary file, reads it back
%   with geometry.receiver and deletes it (nothing is written to the
%   current folder). The segment starts at (x0, z0), dips at dip and is
%   split into Nfault patches of equal width.
%
% INPUTS
%   earthModel   geometry.LDhs object: elastic half-space
%   x0           horizontal position of the up-dip end (m)
%   z0           depth of the up-dip end, positive downward (m)
%   dip          dip angle (degrees)
%   Fault_width  down-dip length of the fault (m)
%   Nfault       number of patches
%
% OUTPUTS
%   rcv          geometry.receiver object: the megathrust patches
%
% NOTES
%   - The file is written with eight columns (n, Vpl = 1, x, z, Width, Dip,
%     W0, qW = 1) and no boundary-condition columns.
%
% SEE ALSO
%   geometry.receiver
%
% AUTHORS
%   Rishav Mallick, 2023

% one-segment fault file, in a temporary location
patchfname = [tempname '.seg'];
removeFile = onCleanup(@() deleteIfPresent(patchfname));

% width of patch segments
w = Fault_width/Nfault;

fileID = fopen(patchfname,'w');
fprintf(fileID,'%s\n',...
    '# n  Vpl,  x   z    Width   Dip   W0   qW');
fprintf(fileID,'%d %.18f %.18f %.18f %.18f %.18f %.18f %d\n',...
    1, 1, x0, z0, Fault_width, dip, w, 1.0);
fclose(fileID);

rcv = geometry.receiver(patchfname,earthModel);

end

function deleteIfPresent(f)
if isfile(f)
    delete(f);
end
end
