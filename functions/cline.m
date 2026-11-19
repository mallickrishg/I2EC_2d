function h = cline(x, y, z, cdata)
% CLINE  Draw a colour-coded line using the edges of a patch.
%
%   h = cline(x, y)
%   h = cline(x, y, z)
%   h = cline(x, y, z, cdata)
%
%   Draws a 2D or 3D line whose colour varies along it, by drawing a patch
%   with no face colour and interpolated edge colours. Used by plotpatch2d
%   to colour fault patches by value.
%
% INPUTS
%   x       vector of x-values
%   y       vector of y-values
%   z       vector of z-values (optional; [] means z = 0)
%   cdata   vector of colour data (optional)
%
% OUTPUTS
%   h       handle to the line (a patch object), returned only if requested
%
% NOTES
%   - 2 inputs: cdata = y, z = 0. 3 inputs: cdata = z. 4 inputs with
%     z = []: z = 0 and cdata as given.
%   - Marker properties ('marker', 'markersize', 'markeredgecolor',
%     'markerfacecolor') work as for a regular line. Line properties
%     ('linestyle', 'linewidth') have no effect on screen but change the
%     output when printing to file; this comes from how MATLAB draws
%     interpolated patch edges. The edge width is set to 3.
%
% EXAMPLE
%   t = 2*pi:.1:8*pi;
%   cline(sqrt(t).*sin(t), sqrt(t).*cos(t)); view(3)                    % Example 1
%   cline(sqrt(t).*sin(t), sqrt(t).*cos(t), t); view(3)                 % Example 2
%   cline(sqrt(t).*sin(t), sqrt(t).*cos(t), t, rand(size(t))); view(3)  % Example 3
%   cline(sqrt(t).*sin(t), sqrt(t).*cos(t), [], rand(size(t))); view(3) % Example 4
%
% SEE ALSO
%   patch, plotpatch2d
%
% AUTHORS
%   S. Hölz, TU-Berlin, seppel_mit_ppATweb.de
%   V 1.0, 16.4.2007
%   Created using Matlab 7.0.4 (SP2)


    % Check input arguments
    error(nargchk(2, 4, nargin))
    if ~isnumeric(x) || ~isnumeric(y) || ~isvector(x) || ~isvector(y) || length(x)~=length(y); 
        error('x and y must be numeric and conforming vectors'); 
    end
    if (nargin == 3 && (~isnumeric(z) || ~isvector(z) || length(x)~=length(z))) || ...
       (nargin == 4 && ~isempty(z) && (~isnumeric(z) || ~isvector(z) || length(x)~=length(z))) ...
       (nargin == 4 && (~isnumeric(cdata) || ~isvector(cdata) || length(x)~=length(cdata)))
        error('z (and cdata) must be a numeric vector and conforming to x and y'); 
    end
    
    % Draw line as patch
    if nargin == 2
        p = patch([x(:)' nan], [y(:)' nan], 0);
        cdata = [y(:)' nan];
    elseif nargin == 3
        p = patch([x(:)' nan], [y(:)' nan], [z(:)' nan], 0);
        cdata = [z(:)' nan];
    elseif nargin == 4
        if isempty(z); z = zeros(size(x)); end
        p = patch([x(:)' nan], [y(:)' nan], [z(:)' nan], 0);
        cdata = [cdata(:)' nan];
    end
    
    set(p,'cdata', cdata, 'edgecolor','interp','facecolor','none','linewidth',3)
    
    
    % Create output
    if nargout == 1; h = p; end