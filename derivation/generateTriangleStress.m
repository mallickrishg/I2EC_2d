function generateTriangleStress()
% GENERATETRIANGLESTRESS  Write functions/computeStressPlaneStrainTriangleShearZone.m.
%
%   generateTriangleStress()
%
%   Differentiates the closed-form displacement of a strained triangle
%   (the edge integrals I22, I23, I32, I33 of
%   computeDisplacementPlaneStrainTriangleShearZone, every branch) with
%   respect to x2 and x3 with the Symbolic Math Toolbox, and the
%   half-space Green's function used by
%   computeDisplacementPlaneStrainTriangleShearZoneTanhSinh. Writes the
%   derivatives as vectorised local functions after the hand-written part
%   in triangleStressTemplate.m. Takes about a minute.
%
% NOTES
%   - The branch expressions are read from the displacement function by
%     line number (the table at the top of this file). If that file
%     changes, update the table and regenerate.
%   - Every atan(u) is replaced by a placeholder whose derivative is
%     (N'D - ND')/(N^2 + D^2), with [N, D] = numden(u): smooth across
%     D = 0, unlike u'/(1 + u^2).
%
% AUTHORS
%   Rishav Mallick, 2026

here=fileparts(mfilename('fullpath'));
root=fileparts(here);
srcFile=fullfile(root,'functions','computeDisplacementPlaneStrainTriangleShearZone.m');
outFile=fullfile(root,'functions','computeStressPlaneStrainTriangleShearZone.m');

% branch: name, line ranges in the displacement file (v: a2 = 0,
% g: general, h: a3 = 0 for I23 and I32)
branches={ ...
    'I22_v',[148 157]; 'I22_g',[160 200]; ...
    'I33_v',[223 244]; 'I33_g',[248 308]; ...
    'I23_v',[332 340]; 'I23_g',[343 347; 351 356; 361 374]; 'I23_h',[343 347; 358 358; 361 374]; ...
    'I32_v',[396 403]; 'I32_g',[406 407; 411 414; 420 435]; 'I32_h',[406 407; 416 416; 420 435]};

src=splitlines(fileread(srcFile));
tmp=tempname; mkdir(tmp); addpath(tmp);
cleanup=onCleanup(@() removeTemp(tmp));
syms x2 x3 m2 m3 a2 a3 t nu y2 y3 real
global ATANREG %#ok<GVMIS>

out=strtrim(fileread(fullfile(here,'triangleStressTemplate.m')));
for b=1:size(branches,1)
    name=branches{b,1}; R=branches{b,2};
    body={};
    for r=1:size(R,1)
        body=[body; src(R(r,1):R(r,2))]; %#ok<AGROW>
    end
    body=strrep(strjoin(body,newline),'atan(','atanPlaceholder(');
    fn=['body_' name];
    fid=fopen(fullfile(tmp,[fn '.m']),'w');
    fprintf(fid,'function I=%s(x2,x3,m2,m3,a2,a3,t,nu)\n%s\nI=us+ui;\nend\n',fn,body);
    fclose(fid);
    writeAtanPlaceholder(tmp);
    rehash

    ATANREG={};
    I=feval(fn,x2,x3,m2,m3,a2,a3,t,nu);
    d2=diff(I,x2); d3=diff(I,x3);
    for k=1:numel(ATANREG)
        u=ATANREG{k}.u; F=ATANREG{k}.F;
        [N,D]=numden(u);
        g2=(diff(N,x2)*D-N*diff(D,x2))/(N^2+D^2);
        g3=(diff(N,x3)*D-N*diff(D,x3))/(N^2+D^2);
        d2=subs(subs(d2,diff(F,x2),g2),diff(F,x3),g3);
        d3=subs(subs(d3,diff(F,x2),g2),diff(F,x3),g3);
        d2=subs(d2,F,atan(u)); d3=subs(d3,F,atan(u));
    end
    f=fullfile(tmp,['d' name '.m']);
    matlabFunction(d2,d3,'Vars',{x2,x3,m2,m3,a2,a3,t,nu},'Outputs',{'d2','d3'}, ...
        'File',f,'Optimize',true,'Comments',sprintf('x2- and x3-derivatives of %s',name));
    out=[out newline newline cleanGenerated(fileread(f))]; %#ok<AGROW>
    fprintf('%s: %d atan terms\n',name,numel(ATANREG));
end

% Green's function (as in computeDisplacementPlaneStrainTriangleShearZoneTanhSinh;
% atan2 instead of atan of a ratio, same for x3 + y3 >= 0, so that the
% derivative is finite at x3 = y3 = 0)
r1=sqrt((x2-y2)^2+(x3-y3)^2); r2=sqrt((x2-y2)^2+(x3+y3)^2);
G22=-1/(2*pi*(1-nu))*((3-4*nu)/4*log(r1)+(8*nu^2-12*nu+5)/4*log(r2)+(x3-y3)^2/r1^2/4 ...
    +((3-4*nu)*(x3+y3)^2+2*y3*(x3+y3)-2*y3^2)/r2^2/4-(y3*x3*(x3+y3)^2)/r2^4);
G23=1/(2*pi*(1-nu))*((1-2*nu)*(1-nu)*atan2(x2-y2,x3+y3)+((x3-y3)*(x2-y2))/r1^2/4 ...
    +(3-4*nu)*((x2-y2)*(x3-y3))/r2^2/4-(y3*x3*(x2-y2)*(x3+y3))/r2^4);
G32=1/(2*pi*(1-nu))*(-(1-2*nu)*(1-nu)*atan2(x2-y2,x3+y3)+(x3-y3)*(x2-y2)/r1^2/4 ...
    +(3-4*nu)*(x2-y2)*(x3-y3)/r2^2/4+y3*x3*(x2-y2)*(x3+y3)/r2^4);
G33=1/(2*pi*(1-nu))*(-(3-4*nu)/4*log(r1)-(8*nu^2-12*nu+5)/4*log(r2) ...
    -(x2-y2)^2/r1^2/4+(2*y3*x3-(3-4*nu)*(x2-y2)^2)/r2^2/4-y3*x3*(x2-y2)^2/r2^4);
Gs=[G22 G23 G32 G33];
dG=num2cell([diff(Gs,x2) diff(Gs,x3)]);
f=fullfile(tmp,'dGreenParts.m');
matlabFunction(dG{:},'Vars',{x2,x3,y2,y3,nu},'File',f,'Optimize',true, ...
    'Outputs',{'g1','g2','g3','g4','g5','g6','g7','g8'}, ...
    'Comments','x2- then x3-derivatives of G22, G23, G32, G33 (point x, source y)');
out=[out newline newline cleanGenerated(fileread(f)) newline newline ...
    'function g=dGreen(x2,x3,y2,y3,nu)' newline ...
    '% derivatives of the Green''s function, as a cell array (see dGreenParts)' newline ...
    'g=cell(1,8);' newline ...
    '[g{:}]=dGreenParts(x2,x3,y2,y3,nu);' newline 'end' newline];

fid=fopen(outFile,'w'); fprintf(fid,'%s',out); fclose(fid);
fprintf('wrote %s\n',outFile);
end

function removeTemp(tmp)
rmpath(tmp); rmdir(tmp,'s');
end

function s=cleanGenerated(s)
% drop the toolbox's version and date lines so regeneration is reproducible
s=regexprep(s,'%\s+This function was generated[^\n]*\n','');
s=regexprep(s,'%\s+\d{2}-\w{3}-\d{4} [\d:]+\n','');
s=strtrim(s);
if ~endsWith(s,'end'), s=[s newline 'end']; end
end

function writeAtanPlaceholder(dir)
f=fullfile(dir,'atanPlaceholder.m');
if isfile(f), return, end
fid=fopen(f,'w');
fprintf(fid,['function y=atanPlaceholder(u)\n' ...
    '%% atan(u) as a symbolic function F_k(x2,x3), recorded with its argument\n' ...
    'global ATANREG\n' ...
    'syms x2 x3 real\n' ...
    'if ~has(u,x2) && ~has(u,x3), y=atan(u); return, end\n' ...
    'k=numel(ATANREG)+1;\n' ...
    'F=str2sym(sprintf(''ATANF%%d(x2,x3)'',k));\n' ...
    'ATANREG{k}=struct(''u'',u,''F'',F);\n' ...
    'y=F;\n' ...
    'end\n']);
fclose(fid);
end
