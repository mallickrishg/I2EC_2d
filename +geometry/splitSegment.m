function [flt]=splitSegment(xo,W,dip,wo,alphaw)
% SPLITSEGMENT  Split a straight segment into patches of growing width.
%
%   flt = geometry.splitSegment(xo, W, dip, wo, alphaw)
%
%   Divides a segment of down-dip length W, starting at xo, into patches
%   whose widths grow geometrically, wo*alphaw^k for k = 0, 1, 2, and so
%   on. When the next width would exceed half of the remaining length,
%   the remainder becomes the last patch, so the widths always sum to W.
%
% INPUTS
%   xo      [2 x 1] upper end of the segment, [x; z] with z the depth,
%           positive downward (m)
%   W       [1 x 1] down-dip length of the segment (m)
%   dip     [1 x 1] dip (degrees)
%   wo      [1 x 1] approximate width of the first patch (m)
%   alphaw  [1 x 1] geometric growth factor of the widths (1 = uniform)
%
% OUTPUTS
%   flt     [Np x 4] one row per patch: x, z, width, dip, with (x, z) the
%           upper end of the patch (z positive downward, m)
%
% SEE ALSO
%   geometry.fpatch.seg2flt, geometry.receiver
%
% AUTHORS
%   Rishav Mallick, 2023

% create wi
Wc=W;
k=0;
w=0;
while Wc>0
    Wt=wo*alphaw^k;
    if Wt > Wc/2
        Wt = Wc;
    end
    wn=min([Wt,Wc]);
    w=[w; wn];
    k=k+1;
    Wc=Wc-wn;
end
Nw=k;

% strike and dip direction normal vectors
Dv=[cosd(dip); sind(dip)];

flt=[];

% loop in dip direction
k=0;
for j=1:Nw
        
    x=xo + sum(w(1:j))*Dv;
    k=k+1;
    
    flt = [flt; [x',w(j+1),dip]];

end

end




        