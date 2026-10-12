function test_negate_exact_solver
%TEST_NEGATE_EXACT_SOLVER  Verification of the negate + exact solver.
%
% Self-contained tests of negateGrowingModes, egh_modal, analyticSolveIVP
% and computePeriodicInitialConditionBalanced, using I2EC's example mesh
% (inputs/ and kernels/).  Run from the I2EC_2d root:
%   addpath('tests'); test_negate_exact_solver
%
% The forcing is built as in MAIN_analyticsolve_singleeqcycle.m, with the
% negated operator: c = edot_lock - A_neg * y_pl.
%
% Tests:
%   1. IVP vs expm of the augmented forced system (round-off agreement);
%      raw and negated operator as input give the same result
%   2. Periodic: residual against expm + agreement with spin-up by doubling
%      (residual also with the raw operator as input)
%   3. Continuity: no jumps in mantle strain rate and integrated strain
%      along a viscosity scan (relative curvature from second differences)
%   4. Series threshold: modal integrals unchanged when the switch moves
%   5. Consistency: IVP from periodic y0 reproduces y(T) + b = y0, with
%      the negated and the raw operator as input
%   6. Integrator mode (lam = 0): egh_modal limits, IVP and periodic state
%   7. Negation on a small hand-built operator with a known answer
%   8. Steady state: with long-term loading only, y_pl stays y_pl
%   9. Balancing: same results as an unbalanced eigendecomposition
%
% No dependency outside I2EC_2d; no absolute paths.
%
% Author: Sharadha Sathiakumar, EOS, 2026

fprintf('\n=== test_negate_exact_solver ===\n');
npass = 0; nfail = 0;

%% ---- Setup: build operator from I2EC example ----
here = fileparts(mfilename('fullpath'));
root = fullfile(here, '..');
addpath(fullfile(root, 'functions'));
addpath(root);   % for +geometry package

yr = 3.15e7;     % seconds per year, as in I2EC's MAIN scripts
nu = 0.25; mu = 30e3;
Trecur = 500 * yr;  Vpl = 1e-9;

earthModel = geometry.LDhs(mu, nu);
rcv = geometry.receiver(fullfile(root, 'inputs', 'megathrust2d.seg'), earthModel);
shz = geometry.shearZoneReceiver(fullfile(root, 'inputs', 'shearzone'), earthModel);
load(fullfile(root, 'kernels', 'evl_orig.mat'), 'evl_orig');

% Rheology (I2EC example)
rcv.Asigma = 1e-6 .* (5e17 / (sum(~rcv.pinnedPosition .* rcv.W))) .* ones(rcv.N, 1);
oceanic = shz.xc(:,1) < -shz.xc(:,2) / tand(rcv.dip(1));
eta_mantle = 1e19 * ones(shz.N, 1);
eta_mantle(~oceanic) = 10^18.5;
locked = abs(rcv.xc(:,2)) > 0 & abs(rcv.xc(:,2)) < 30e3;
rcv.Vpl = Vpl .* ones(rcv.N, 1);
[e22_dev_lt, e23_lt] = getStrainratesLongterm(shz, rcv.dip(1)*pi/180, [0, 20e3], [-140e3, 35e3], 1);
y_pl = [rcv.Vpl(~locked); e22_dev_lt .* Vpl; -e23_lt .* Vpl];

stresskernel = assembleStressKernel(evl_orig, locked);
slip_co = zeros(rcv.N, 1);
slip_co(locked) = Trecur * Vpl;
flv = zeros(rcv.N, 1);  flv(locked) = Vpl;
sys = @(fac) build_system(evl_orig, stresskernel, locked, rcv.Asigma, eta_mantle * fac, slip_co, flv, y_pl);

[A, A_neg, c, b, c_lock] = sys(1);
n = size(A, 1);
imant = (n - 2*shz.N + 1):n;    % mantle strain-rate components

lam0 = eig(A);
fprintf('System size n = %d, growing modes negated: %d\n\n', n, sum(real(lam0) > 0));

y0_per = computePeriodicInitialConditionBalanced(A_neg, c, b, Trecur);
tvec = [40/365, 1, 10, 50, 500] * yr;     % 40 d to 500 yr

%% ---- Test 1: IVP vs expm of augmented system ----
fprintf('Test 1: IVP vs expm ... ');

[Y, U] = analyticSolveIVP(A_neg, c, y0_per, tvec);
[Yref, Uref] = expm_ivp(A_neg, c, y0_per, tvec);
err_y = max(vecnorm(Y - Yref) ./ vecnorm(Yref));
err_u = max(vecnorm(U - Uref) ./ vecnorm(Uref));

% raw operator as input: the solver negates it itself
[Yr, Ur] = analyticSolveIVP(A, c, y0_per, tvec);
err_raw = max([vecnorm(Yr - Y) ./ vecnorm(Y), vecnorm(Ur - U) ./ vecnorm(U)]);

tol1 = 1e-10;
pass1 = max([err_y, err_u, err_raw]) < tol1;
report(pass1, sprintf('max rel err: y %.2e, u %.2e; raw vs negated input %.2e (tol %.0e)', ...
    err_y, err_u, err_raw, tol1));
[npass, nfail] = tally(pass1, npass, nfail);

%% ---- Test 2: Periodic residual + spin-up ----
fprintf('Test 2: periodic residual + spin-up ... ');

% 2a: residual of (I - Phi) y0 = eta + b, Phi and eta from expm
eM = expm([A_neg c; zeros(1, n+1)] * Trecur);
Phi = eM(1:n, 1:n);
eta = eM(1:n, n+1);
res_rel = norm((eye(n) - Phi) * y0_per - (eta + b)) / norm(eta + b);

% raw operator as input: the solver negates it itself
y0_raw = computePeriodicInitialConditionBalanced(A, c, b, Trecur);
res_rel = max(res_rel, norm((eye(n) - Phi) * y0_raw - (eta + b)) / norm(eta + b));

% 2b: spin-up by doubling (2^40 cycles from y = 0)
F = Phi;  s = eta + b;
for k = 1:40
    s = F * s + s;
    F = F * F;
end
spinup_rel = norm(y0_per - s) / norm(y0_per);
rho = max(abs(eig(Phi)));

tol2a = 1e-10;
tol2b = 1e-8;   % spin-up accumulates round-off over 40 matrix squarings
pass2 = res_rel < tol2a && spinup_rel < tol2b && rho < 1;
report(pass2, sprintf('residual %.2e (tol %.0e), spin-up %.2e (tol %.0e), rho(Phi) %.4f', ...
    res_rel, tol2a, spinup_rel, tol2b, rho));
[npass, nfail] = tally(pass2, npass, nfail);

%% ---- Test 3: Continuity along a viscosity scan ----
fprintf('Test 3: continuity scan ... ');

% Mantle viscosity times 10^(-0.5 ... 0.5).  Over this range eigenvalues
% of the raw operator cross the old solver's thresholds (|Re| = 1e-10),
% where the old solver jumps.  Quantities: mantle strain rate y and
% integrated strain u (periodic cycle, at 10 and 100 yr).  For each, the
% central second difference along the scan, relative to the quantity's
% size and divided by the step^2 (log10 units), estimates the relative
% curvature, which does not depend on the number of scan points.
npts = 41;
log_factors = linspace(-0.5, 0.5, npts);
t3 = [10, 100] * yr;
Q = zeros(4 * 2*shz.N, npts);
for ip = 1:npts
    [~, An, cn, bn] = sys(10^log_factors(ip));
    y0n = computePeriodicInitialConditionBalanced(An, cn, bn, Trecur);
    [Yn, Un] = analyticSolveIVP(An, cn, y0n, t3);
    Q(:, ip) = [Yn(imant, 1); Yn(imant, 2); Un(imant, 1); Un(imant, 2)];
end
nq = 2*shz.N;
ds = log_factors(2) - log_factors(1);
curv = 0;
for iq = 1:4
    q = Q((iq-1)*nq + (1:nq), :);
    d2 = vecnorm(q(:, 3:end) - 2*q(:, 2:end-1) + q(:, 1:end-2)) ./ vecnorm(q(:, 2:end-1));
    curv = max(curv, max(d2) / ds^2);
end

% The smooth solution has a relative curvature of about 5 per decade^2;
% a jump of relative size J gives about J/ds^2 (the old solver: ~350).
tol3 = 10;
pass3 = curv < tol3;
report(pass3, sprintf('max relative curvature %.2f per decade^2 (tol %g)', curv, tol3));
[npass, nfail] = tally(pass3, npass, nfail);

%% ---- Test 4: Series threshold sensitivity ----
fprintf('Test 4: series threshold ... ');

lam_neg = eig(A_neg);
max_diff = 0;
for t4 = [40/365, 1, 10, 100] * yr
    [~, g_def, h_def] = egh_modal(lam_neg, t4);
    for thr = [1e-3, 1e-1]
        [~, g_t, h_t] = egh_thresh(lam_neg, t4, thr);
        max_diff = max([max_diff, max(abs(g_t - g_def)) / max(abs(g_def)), ...
                        max(abs(h_t - h_def)) / max(abs(h_def))]);
    end
end

tol4 = 1e-11;   % at a switch of 1e-3, closed-form h is accurate to ~eps/|x| = 2e-13
pass4 = max_diff < tol4;
report(pass4, sprintf('max rel diff in g, h %.2e (tol %.0e)', max_diff, tol4));
[npass, nfail] = tally(pass4, npass, nfail);

%% ---- Test 5: IVP-periodic consistency ----
fprintf('Test 5: IVP-periodic consistency ... ');

Y_T = analyticSolveIVP(A_neg, c, y0_per, Trecur);
cycle_err = norm(Y_T + b - y0_per) / norm(y0_per);

% raw operator as input to both solvers
Y_Tr = analyticSolveIVP(A, c, y0_raw, Trecur);
cycle_err = max(cycle_err, norm(Y_Tr + b - y0_raw) / norm(y0_raw));

tol5 = 1e-10;
pass5 = cycle_err < tol5;
report(pass5, sprintf('|y(T) + b - y0| / |y0| = %.2e (tol %.0e)', cycle_err, tol5));
[npass, nfail] = tally(pass5, npass, nfail);

%% ---- Test 6: integrator mode (lam = 0) ----
fprintf('Test 6: integrator mode (lam = 0) ... ');

% 6a: egh_modal at lam = 0 and |lam| far below the series switch
t6 = 10 * yr;
[e6, g6, h6] = egh_modal([0; 1e-30; -1e-30i], t6);
err6a = max([abs(e6 - 1); abs(g6 - t6)/t6; abs(h6 - t6^2/2)/(t6^2/2)]);

% 6b: IVP with an exact integrator appended, dy/dt = c0 for the extra
% component: y = y00 + c0 t, u = y00 t + c0 t^2/2; the other components
% must equal the solution without it
c0 = mean(abs(c));  y00 = mean(abs(y0_per));
A6 = blkdiag(A_neg, 0);
[Y6, U6] = analyticSolveIVP(A6, [c; c0], [y0_per; y00], tvec);
y_int = y00 + c0 * tvec;
u_int = y00 * tvec + c0 * tvec.^2 / 2;
err6b = max([abs(Y6(end,:) - y_int) ./ abs(y_int), abs(U6(end,:) - u_int) ./ abs(u_int)]);
err6b_rest = max([vecnorm(Y6(1:n,:) - Y) ./ vecnorm(Y), vecnorm(U6(1:n,:) - U) ./ vecnorm(U)]);

% 6c: periodic with the integrator; a periodic solution exists only if
% c0*T + b0 = 0, and the guard |lam| < 1e-16 gives the minimum-norm z0 = 0
y0_6 = computePeriodicInitialConditionBalanced(A6, [c; c0], [b; -c0*Trecur], Trecur);
err6c = abs(y0_6(end)) / norm(y0_per);
err6c_rest = norm(y0_6(1:n) - y0_per) / norm(y0_per);

tol6 = 1e-10;
pass6 = max([err6a, err6b, err6b_rest, err6c, err6c_rest]) < tol6;
report(pass6, sprintf(['egh %.1e; IVP integrator %.1e, rest %.1e; ' ...
    'periodic integrator %.1e, rest %.1e (tol %.0e)'], ...
    err6a, err6b, err6b_rest, err6c, err6c_rest, tol6));
[npass, nfail] = tally(pass6, npass, nfail);

%% ---- Test 7: negation on a hand-built operator ----
fprintf('Test 7: hand-built negation ... ');

% A = S D S^-1 with a real S and a block-diagonal real D holding a stable
% mode, a slow stable mode, a growing real mode and a growing complex pair
% (block B = [a w; -w a], eigenvalues a +- i w).  Negating every eigenvalue
% with Re > 0 while keeping the eigenvectors turns B into -B, so the
% expected operator is S D_neg S^-1, built without any eigendecomposition.
rng(0);
S7 = eye(5) + 0.3 * randn(5);
a = 2e-10;  w = 5e-10;  B = [a w; -w a];
D7     = blkdiag(-1e-9, -1e-13, 3e-10, B);
D7_neg = blkdiag(-1e-9, -1e-13, -3e-10, -B);
A7 = S7 * D7 / S7;
A7_ref = S7 * D7_neg / S7;
c7 = randn(5, 1) * 1e-10;  y7 = randn(5, 1);  b7 = randn(5, 1);
t7 = [1e6, 1e8, 1e9, 1e10];     % |lam t| from 1e-7 (series) to 10

[~, ~, A7_neg] = negateGrowingModes(A7);
err7a = norm(A7_neg - A7_ref) / norm(A7_ref);
[Y7, U7] = analyticSolveIVP(A7, c7, y7, t7);
[Y7r, U7r] = expm_ivp(A7_ref, c7, y7, t7);
err7b = max([vecnorm(Y7 - Y7r) ./ vecnorm(Y7r), vecnorm(U7 - U7r) ./ vecnorm(U7r)]);
T7 = 1e10;
y07 = computePeriodicInitialConditionBalanced(A7, c7, b7, T7);
eM7 = expm([A7_ref c7; zeros(1, 6)] * T7);
err7c = norm((eye(5) - eM7(1:5, 1:5)) * y07 - (eM7(1:5, 6) + b7)) / norm(eM7(1:5, 6) + b7);

tol7 = 1e-10;
pass7 = max([err7a, err7b, err7c]) < tol7;
report(pass7, sprintf('operator %.1e, IVP vs expm %.1e, periodic residual %.1e (tol %.0e)', ...
    err7a, err7b, err7c, tol7));
[npass, nfail] = tally(pass7, npass, nfail);

%% ---- Test 8: long-term state is the steady state ----
fprintf('Test 8: steady state ... ');

% With long-term loading only (no locking, no earthquake), the system
% started at y_pl must stay there: y = y_pl, u = y_pl t; the periodic
% state with b = 0 is y_pl.  This holds because c is built with the
% negated operator.
c8 = c - c_lock;                 % = -A_neg * y_pl
[Y8, U8] = analyticSolveIVP(A_neg, c8, y_pl, tvec);
err8a = max([vecnorm(Y8 - y_pl) / norm(y_pl), vecnorm(U8 - y_pl * tvec) ./ (norm(y_pl) * tvec)]);
y08 = computePeriodicInitialConditionBalanced(A_neg, c8, zeros(n, 1), Trecur);
err8b = norm(y08 - y_pl) / norm(y_pl);

% for information: forcing built with the raw operator instead
Y8raw = analyticSolveIVP(A, -A * y_pl, y_pl, tvec(end));
drift_raw = norm(Y8raw(imant) - y_pl(imant)) / norm(y_pl(imant));

tol8 = 1e-10;
pass8 = max(err8a, err8b) < tol8;
report(pass8, sprintf(['IVP %.1e, periodic %.1e (tol %.0e); ' ...
    'with raw-operator forcing the mantle strain rate drifts by %.1f%% in 500 yr'], ...
    err8a, err8b, tol8, 100 * drift_raw));
[npass, nfail] = tally(pass8, npass, nfail);

%% ---- Test 9: balanced vs unbalanced eigendecomposition ----
fprintf('Test 9: balanced vs unbalanced ... ');

[Vu, Du] = eig(A);
lu = diag(Du);
lu(real(lu) > 0) = -lu(real(lu) > 0);
z0u = Vu \ y0_per;  ctu = Vu \ c;  zbu = Vu \ b;
Yu = zeros(n, numel(tvec));  Uu = Yu;
for j = 1:numel(tvec)
    [ej, gj, hj] = egh_modal(lu, tvec(j));
    Yu(:,j) = real(Vu * (ej .* z0u + gj .* ctu));
    Uu(:,j) = real(Vu * (gj .* z0u + hj .* ctu));
end
[~, gT] = egh_modal(lu, Trecur);
y0u = real(Vu * ((gT .* ctu + zbu) ./ (-lu .* gT)));
[Yb, Ub] = analyticSolveIVP(A, c, y0_per, tvec);
y0b = computePeriodicInitialConditionBalanced(A, c, b, Trecur);
err9 = max([vecnorm(Yb - Yu) ./ vecnorm(Yu), vecnorm(Ub - Uu) ./ vecnorm(Uu), ...
            norm(y0b - y0u) / norm(y0u)]);

tol9 = 1e-9;
pass9 = err9 < tol9;
report(pass9, sprintf('max rel diff (y, u, periodic y0) %.2e (tol %.0e)', err9, tol9));
[npass, nfail] = tally(pass9, npass, nfail);

%% ---- Summary ----
fprintf('\n=== %d passed, %d failed ===\n\n', npass, nfail);
if nfail > 0
    error('test_negate_exact_solver: %d test(s) FAILED', nfail);
end

end

%% ======== Local functions ========

function [A, A_neg, c, b, c_lock] = build_system(evl_orig, stresskernel, locked, Asigma, eta_mantle, slip_co, flv, y_pl)
% Operator, forcing and coseismic jump as in MAIN_analyticsolve_singleeqcycle.m
viscvec = [Asigma(~locked); eta_mantle * 1e-6; eta_mantle * 1e-6];
A = diag(1 ./ viscvec) * stresskernel;
[~, ~, A_neg] = negateGrowingModes(A);
b = [evl_orig.KK(~locked,:) * slip_co; evl_orig.KL(:,:,1) * slip_co; evl_orig.KL(:,:,2) * slip_co] ./ viscvec;
c_lock = -[evl_orig.KK(~locked,:) * flv; evl_orig.KL(:,:,1) * flv; evl_orig.KL(:,:,2) * flv] ./ viscvec;
c = c_lock - A_neg * y_pl;
end

function [Y, U] = expm_ivp(A, c, y0, tvec)
% Reference: expm of the augmented system d/dt [u; y; 1] = [0 I 0; 0 A c; 0 0 0] [u; y; 1]
n = size(A, 1);
M = [zeros(n) eye(n) zeros(n,1); zeros(n) A c; zeros(1, 2*n+1)];
Y = zeros(n, numel(tvec));  U = Y;
for j = 1:numel(tvec)
    r = expm(M * tvec(j)) * [zeros(n,1); y0; 1];
    U(:,j) = r(1:n);
    Y(:,j) = r(n+1:2*n);
end
end

function report(pass, msg)
if pass, fprintf('PASS  %s\n', msg);
else,    fprintf('FAIL  %s\n', msg);
end
end

function [npass, nfail] = tally(pass, npass, nfail)
if pass, npass = npass + 1; else, nfail = nfail + 1; end
end

function [e, g, h] = egh_thresh(lam, t, thresh)
% egh_modal with a configurable series switch (for test 4 only)
x = lam * t; e = exp(x); g = zeros(size(lam)); h = g;
s = abs(x) < thresh;
xs = x(s); gs = zeros(size(xs)); hs = gs; fk = 1;
for k = 0:7
    gs = gs + xs.^k / fk;
    hs = hs + xs.^k / (fk * (k+2));
    fk = fk * (k+2);
end
g(s) = t * gs; h(s) = t^2 * hs;
xb = x(~s); em1 = expm1(xb);
g(~s) = em1 ./ lam(~s);
h(~s) = (em1 - xb) ./ lam(~s).^2;
end
