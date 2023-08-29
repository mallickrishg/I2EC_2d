function rcv = create_megathrust(earthModel,x0,z0,dip,Fault_width,Nfault)

% megathrust fault file
patchfname = '../megathrust2d.seg';

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