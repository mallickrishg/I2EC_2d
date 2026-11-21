% COUPLEDEXAMPLE  Fault-mantle gap and coupled-operator check on the example model.
%
% Measures how far each megathrust patch is from the mantle mesh, checks
% Betti reciprocity between KL and LK (W*LK = 2*a*KL'), and certifies the
% width/area-weighted coupled operator [KK LK; KL LL]. Writes
% fault_mantle_gap.png next to this script. See tasks.md, section 0a.
here = fileparts(mfilename('fullpath'));
root = fileparts(fileparts(here));
addpath(fullfile(root, 'examples'));
m = exampleModel('StressKernels', true);
rcv = m.rcv; shz = m.shz; evl = m.evl;
N = rcv.N; M = shz.N;

%% geometry: fault patch endpoints
xc = rcv.xc; W = rcv.W(:); dv = rcv.dv;          % dv = down-dip unit vector?
fprintf('rcv.xc z range: %.1f .. %.1f km; shz.xc z range %.1f .. %.1f km\n', ...
    min(xc(:,2))/1e3, max(xc(:,2))/1e3, min(shz.xc(:,2))/1e3, max(shz.xc(:,2))/1e3);
A = xc - 0.5*W.*dv; B = xc + 0.5*W.*dv;          % patch ends
% distance from points on the fault (dense) to the mantle mesh (triangles)
t = linspace(0, 1, 50);
P = []; owner = [];
for i = 1:N
    P = [P; A(i,:) + t'.*(B(i,:) - A(i,:))]; owner = [owner; i*ones(numel(t),1)]; %#ok<AGROW>
end
tri = shz.tri; V = shz.vert;
dmin = inf(size(P,1),1);
for k = 1:size(tri,1)
    T = V(tri(k,:),:);
    dmin = min(dmin, pointTriangleDistance(P, T));
end
dpatch = accumarray(owner, dmin, [N 1], @min);
dcent = inf(N,1);
for k = 1:size(tri,1)
    dcent = min(dcent, pointTriangleDistance(xc, V(tri(k,:),:)));
end
fprintf('fault-to-mantle distance: min over fault %.2f km (patch %d), median over patches %.1f km\n', ...
    min(dpatch)/1e3, find(dpatch == min(dpatch), 1), median(dpatch)/1e3);
fprintf('  per patch (km): %s\n', mat2str(round(dpatch'/1e3, 1)));
fprintf('  collocation-point distance min %.2f km\n', min(dcent)/1e3);
% size of nearest cells relative to distance
h = sqrt(2*shz.area(:));
[~, near] = min(vecnorm(shz.xc - xc(find(dpatch == min(dpatch),1),:), 2, 2));
fprintf('nearest cell size ~%.1f km; fault patch width %.1f km\n', h(near)/1e3, W(1)/1e3);

%% Betti weighting: W_i LK(i,j,c) vs alpha * a_j KL(j,i,c)
a = shz.area(:);
for c = 1:2
    X = W .* evl.LK(:,:,c);            % N x M
    Y = (a .* evl.KL(:,:,c))';         % N x M
    alpha = (X(:)'*Y(:)) / (Y(:)'*Y(:));
    fprintf('component %d: W*LK ~ %.4f * (a*KL)''  ; rel. misfit %.3e\n', c, alpha, ...
        norm(X - alpha*Y, 'fro')/norm(X, 'fro'));
end
% KK symmetry and definiteness (width-weighted)
SK = W .* evl.KK;
fprintf('KK: width-weighted asymmetry %.3e, max eig of sym part / max|eig| = %.3e\n', ...
    norm(SK - SK', 'fro')/norm(SK, 'fro'), max(eig((SK+SK')/2))/max(abs(eig((SK+SK')/2))));

%% coupled operator, weighted by work: fault rows by W, mantle rows by 2a
K = assembleStressKernel(evl, false(N,1));
S = [W; 2*a; 2*a] .* K;
dd = 1 ./ sqrt(abs(diag(S)));            % Jacobi scaling (congruence keeps signs)
St = dd .* S .* dd';
lam = eig((St + St')/2);
C1 = St(1:N, N+1:end); C2 = St(N+1:end, 1:N)';
mis = vecnorm(C1 - C2, 2, 2) ./ vecnorm(C1, 2, 2);
fprintf('coupled operator (scaled): sym part max eig %.3e, %d >= 0; fault block alone max eig %.3e\n', ...
    max(lam), sum(lam >= 0), max(eig((St(1:N,1:N) + St(1:N,1:N)')/2)));
fprintf('KL/LK mismatch per patch: median %.3f; patches %d-%d: %s\n', median(mis), N-3, N, ...
    mat2str(round(mis(end-3:end)', 3)));

%% eigenvalues of K*D for random positive D
rng(1); worst = -inf;
for trial = 1:20
    D = exp(4*randn(N+2*M, 1));
    mu = eig(K .* D');    % K*diag(D)
    worst = max(worst, max(real(mu))/max(abs(mu)));
end
fprintf('K*diag(D), 20 random D (log-normal, 4 e-folds): max Re/max|.| = %.3e\n', worst);


%% figure: fault and mantle mesh
G = struct('A', A, 'B', B);
f = figure('Visible', 'off', 'Position', [0 0 1300 520]);
for p = 1:2
    ax = subplot(1, 2, p); hold on; box on; axis equal
    patch('Faces', shz.tri, 'Vertices', shz.vert/1e3, 'FaceColor', [0.85 0.9 1], 'EdgeColor', [0.45 0.55 0.8]);
    for i = 1:rcv.N
        plot([G.A(i,1) G.B(i,1)]/1e3, [G.A(i,2) G.B(i,2)]/1e3, '-', 'Color', [0.8 0.1 0.1], 'LineWidth', 2);
    end
    plot(rcv.xc(:,1)/1e3, rcv.xc(:,2)/1e3, 'k.', 'MarkerSize', 8);
    xlabel('x (km)'); ylabel('z (km)');
    if p == 1
        xlim([-60 300]); ylim([-120 5]);
        title('Megathrust (red) and mantle mesh: 19.4 km gap below patches 1-45');
    else
        xlim([170 215]); ylim([-62 -38]);
        title('Down-dip end: patches 46-49 are 17, 13, 9, 1 km from the mantle');
        for i = rcv.N-3:rcv.N
            text(rcv.xc(i,1)/1e3 - 2.5, rcv.xc(i,2)/1e3 + 1.2, sprintf('%d', i), 'FontSize', 10);
        end
    end
end
exportgraphics(f, fullfile(here, 'fault_mantle_gap.png'), 'Resolution', 150);

function d = pointTriangleDistance(P, T)
% distance from points P [n x 2] to triangle T [3 x 2] (0 inside)
in = inpolygon(P(:,1), P(:,2), T(:,1), T(:,2));
d = inf(size(P,1),1);
for e = 1:3
    a = T(e,:); b = T(mod(e,3)+1,:);
    ab = b - a; s = ((P - a)*ab') / (ab*ab'); s = min(max(s, 0), 1);
    d = min(d, vecnorm(P - (a + s.*ab), 2, 2));
end
d(in) = 0;
end
