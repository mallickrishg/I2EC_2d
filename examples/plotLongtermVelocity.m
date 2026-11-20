% Plot the long-term (plate-rate) velocity field of the example model:
% the corner-flow solution that loads the earthquake-cycle examples
% (computeLongTermVelocity), on a grid over the model domain.
% Run from any folder:  run examples/plotLongtermVelocity

clear
addpath(fileparts(mfilename('fullpath')))
m = exampleModel();
rcv = m.rcv; boundary = m.boundary; shz = m.shz;

Vpl = m.Vpl;        % m/s
Tplate = 20e3;      % oceanic plate thickness (m)

figure(1),clf
plotpatch2d(rcv, rcv.xc(:,2)./1e3), hold on
plotpatch2d(boundary, boundary.xc(:,2)./1e3)
plotshz2d(shz)
axis tight equal
box on
set(gca, 'YDir', 'normal', 'Fontsize', 20, 'Linewidth', 2)

%% observation points
x = linspace(min(boundary.x(:,1)), max(boundary.x(:,1)), 50);
y = linspace(max(boundary.x(:,2)), min(boundary.x(:,2)), 50);
[X, Y] = meshgrid(x, y);
X = X(:); Y = Y(:); dip = rcv.dip(1);
ocean_translate_params = [0, 20e3];
continental_translate_params = [-140e3, 35e3];

[Vx, Vz] = computeLongTermVelocity(X, Y, Vpl, Tplate, dip, shz, boundary, ...
    ocean_translate_params, continental_translate_params);

hold on
axis equal
plot(X/1e3, Y/1e3, 'c.')
quiver(X/1e3, Y/1e3, Vx, Vz, 'k')
