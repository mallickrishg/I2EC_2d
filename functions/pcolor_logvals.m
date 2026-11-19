function pcolor_logvals(x,y,toplot,minval)
% PCOLOR_LOGVALS  Prepare log-scaled values of signed data (unfinished).
%
%   pcolor_logvals(x, y, toplot, minval)
%
%   Intended to plot positive and negative values with log-scale colours.
%   It clips |toplot| below minval to minval and takes log10 of the
%   positive and negative parts separately.
%
% INPUTS
%   x, y    plotting coordinates (not used yet)
%   toplot  values to colour
%   minval  smallest magnitude kept; smaller magnitudes are set to it
%
% OUTPUTS
%   none
%
% NOTES
%   - Unfinished: the function draws nothing and returns nothing; the log
%     values are computed and discarded. log10 of the negative part is
%     complex. It is not called anywhere in I2EC.
%
% AUTHORS
%   Rishav Mallick, 2023

posindex = toplot>=0;
negindex = ~posindex;

toplot(posindex & toplot < minval) = minval;
toplot(negindex & toplot > -minval) = -minval;

pos_logvals = zeros(size(toplot));
pos_logvals(posindex) = log10(toplot(posindex));

neg_logvals = zeros(size(toplot));
neg_logvals(negindex) = log10(toplot(negindex));




end