function newmap = bluewhitered(m, varargin)
% BLUEWHITERED  Blue, white and red colour map centred on zero.
%
%   newmap = bluewhitered(m)
%   newmap = bluewhitered(m, lims)
%   newmap = bluewhitered
%
%   Returns an m-by-3 colour map that runs from blue through white to red,
%   with white at the colour-axis value closest to zero. Useful for
%   images and surfaces with positive and negative values. If the colour
%   limits are all positive (or all negative), only the white-to-red (or
%   blue-to-white) half is used.
%
% INPUTS
%   m       number of colours (default: length of the current colormap)
%   lims    [1 x 2] colour limits (default: CLim of the current axes)
%
% OUTPUTS
%   newmap  [m x 3] RGB colour map
%
% EXAMPLE
%   figure
%   imagesc(peaks(250));
%   colormap(bluewhitered(256)), colorbar
%
%   figure
%   imagesc(peaks(250), [0 8])
%   colormap(bluewhitered), colorbar
%
%   figure
%   imagesc(peaks(250), [-6 0])
%   colormap(bluewhitered), colorbar
%
%   figure
%   surf(peaks)
%   colormap(bluewhitered)
%   axis tight
%
% SEE ALSO
%   colormap, rgbplot, hsv, hot, cool, bone, copper, pink, flag
%
% AUTHORS
%   Third-party colour map, original author not named in this file,
%   added to I2EC by Rishav Mallick, 2023


if nargin < 1
   m = size(get(gcf,'colormap'),1);
end


bottom = [0 0 0.5];
botmiddle = [0 0.5 1];
middle = [1 1 1];
topmiddle = [1 0 0];
top = [0.5 0 0];

% Find middle
if nargin == 2;
	lims = varargin{:};
else
	lims = get(gca, 'CLim');
end

% Find ratio of negative to positive
if (lims(1) < 0) & (lims(2) > 0)
    % It has both negative and positive
    % Find ratio of negative to positive
    ratio = abs(lims(1)) / (abs(lims(1)) + lims(2));
    neglen = round(m*ratio);
    poslen = m - neglen;
    
    % Just negative
    new = [bottom; botmiddle; middle];
    len = length(new);
    oldsteps = linspace(0, 1, len);
    newsteps = linspace(0, 1, neglen);
    newmap1 = zeros(neglen, 3);
    
    for i=1:3
        % Interpolate over RGB spaces of colormap
        newmap1(:,i) = min(max(interp1(oldsteps, new(:,i), newsteps)', 0), 1);
    end
    
    % Just positive
    new = [middle; topmiddle; top];
    len = length(new);
    oldsteps = linspace(0, 1, len);
    newsteps = linspace(0, 1, poslen);
    newmap = zeros(poslen, 3);
    
    for i=1:3
        % Interpolate over RGB spaces of colormap
        newmap(:,i) = min(max(interp1(oldsteps, new(:,i), newsteps)', 0), 1);
    end
    
    % And put 'em together
    newmap = [newmap1; newmap];
    
elseif lims(1) >= 0
    % Just positive
    new = [middle; topmiddle; top];
    len = length(new);
    oldsteps = linspace(0, 1, len);
    newsteps = linspace(0, 1, m);
    newmap = zeros(m, 3);
    
    for i=1:3
        % Interpolate over RGB spaces of colormap
        newmap(:,i) = min(max(interp1(oldsteps, new(:,i), newsteps)', 0), 1);
    end
    
else
    % Just negative
    new = [bottom; botmiddle; middle];
    len = length(new);
    oldsteps = linspace(0, 1, len);
    newsteps = linspace(0, 1, m);
    newmap = zeros(m, 3);
    
    for i=1:3
        % Interpolate over RGB spaces of colormap
        newmap(:,i) = min(max(interp1(oldsteps, new(:,i), newsteps)', 0), 1);
    end
    
end
% 
% m = 64;
% new = [bottom; botmiddle; middle; topmiddle; top];
% % x = 1:m;
% 
% oldsteps = linspace(0, 1, 5);
% newsteps = linspace(0, 1, m);
% newmap = zeros(m, 3);
% 
% for i=1:3
%     % Interpolate over RGB spaces of colormap
%     newmap(:,i) = min(max(interp1(oldsteps, new(:,i), newsteps)', 0), 1);
% end
% 
% % set(gcf, 'colormap', newmap), colorbar