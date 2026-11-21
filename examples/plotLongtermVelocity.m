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

%% observation points
x = linspace(min(boundary.x(:,1)), max(boundary.x(:,1)), 36);
y = linspace(max(boundary.x(:,2)), min(boundary.x(:,2)), 15);
[X, Y] = meshgrid(x, y);
X = X(:); Y = Y(:); dip = rcv.dip(1);
ocean_translate_params = [0, 20e3];
continental_translate_params = [-140e3, 35e3];

[Vx, Vz] = computeLongTermVelocity(X, Y, Vpl, Tplate, dip, shz, boundary, ...
    ocean_translate_params, continental_translate_params);

%% plot: velocity arrows over the model, mantle mesh shaded
st = exampleStyle();
st.figure(1, 1000, 560)
plotShz2d(shz, ones(shz.N, 1)), hold on
set(findobj(gca, 'Type', 'patch'), 'FaceColor', st.blue(40,:))
plotPatches2d(boundary)
plotPatches2d(rcv)
quiver(X/1e3, Y/1e3, Vx, Vz, 0.9, 'Color', st.ink, 'LineWidth', 0.75, 'MaxHeadSize', 0.25)
xlabel('distance from the trench (km)'), ylabel('depth (km)')
title('Long-term velocity (plate rate, corner flow)')
st.map()
set(findobj(gca, 'Type', 'line'), 'Color', st.ink2, 'LineWidth', 1)   % boundary
plotPatches2d(rcv), set(findobj(gca, 'Type', 'line', '-not', 'Color', st.ink2), 'Color', st.ink, 'LineWidth', 2.5)
