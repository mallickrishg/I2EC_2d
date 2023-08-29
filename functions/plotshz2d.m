function plotshz2d(shz,toplot)
km2m = 1e3;

if nargin == 2
    x = [shz.A(:,1),shz.B(:,1),shz.C(:,1)]';
    y = [shz.A(:,2),shz.B(:,2),shz.C(:,2)]';
    
    patch(x./km2m,y./km2m,toplot)
else
    for i = 1:shz.N
        plot([shz.vert(shz.tri(i,:),1);shz.vert(shz.tri(i,1),1)]./km2m,...
            [shz.vert(shz.tri(i,:),2);shz.vert(shz.tri(i,1),2)]./km2m,'r-'), hold on
    end
end
end