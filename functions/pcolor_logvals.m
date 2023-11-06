function pcolor_logvals(x,y,toplot,minval)
% function to plot positive/negative values with log-scale coloring

posindex = toplot>=0;
negindex = ~posindex;

toplot(posindex & toplot < minval) = minval;
toplot(negindex & toplot > -minval) = -minval;

pos_logvals = zeros(size(toplot));
pos_logvals(posindex) = log10(toplot(posindex));

neg_logvals = zeros(size(toplot));
neg_logvals(negindex) = log10(toplot(negindex));




end