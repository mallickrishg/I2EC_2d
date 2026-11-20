classdef (SharedTestFixtures = {I2ECPathFixture}) GreensFunctionsTest < matlab.unittest.TestCase
    %GREENSFUNCTIONSTEST  Elementary half-space solutions (lowest level of I2EC).
    %
    % Line dislocations (geometry.LDdispHS / geometry.LDstressHS, from
    % CutAndDisplace) and triangular strain volumes (Barbot 2018:
    % computeDisplacementPlaneStrainTriangleShearZone*, computeStress*;
    % the analytic stress against independent derivatives and identities).
    %
    % Every check here is a physical identity that must hold whatever the
    % sign convention: free surface is traction-free, the displacement
    % jumps by the prescribed slip, stress = Hooke(gradient of
    % displacement), and two independent implementations agree.

    properties (Constant)
        G = 30e3
        nu = 0.25
        % line dislocation used throughout: centre, half-length, angle
        LDxe = 3
        LDye = -10
        LDa = 2
        LDbeta = -deg2rad(30)
        % triangle (x2 east, x3 depth positive) and observation points
        TriA = [10 5]
        TriB = [20 12]
        TriC = [8 18]
        ObsX2 = [0 30 15 40 14]   % pt 5 (14,12) is inside the triangle
        ObsX3 = [0  2 25 40 12]
    end

    methods (Test)
        %% ---------------- line dislocation ----------------
        function lineDislocationSurfaceIsTractionFree(tc)
            xs = linspace(-50, 50, 11)';
            E = 2*tc.G*(1 + tc.nu);
            S = geometry.LDstressHS(xs, 0*xs, tc.LDxe, tc.LDye, tc.LDa, tc.LDbeta, 1, 0.5, tc.nu, E);
            scale = max(abs(S(:)));
            tc.verifyLessThan(max(abs(S(:,2))) / scale, 1e-10, 'Syy must vanish at y = 0');
            tc.verifyLessThan(max(abs(S(:,3))) / scale, 1e-10, 'Sxy must vanish at y = 0');
        end

        function lineDislocationJumpEqualsSlip(tc)
            % displacement jump across the element midpoint = prescribed Ds / Dn
            n = [-sin(tc.LDbeta) cos(tc.LDbeta)];      % element normal
            d = 1e-6;
            P = [tc.LDxe tc.LDye] + d*n;
            M = [tc.LDxe tc.LDye] - d*n;
            for DsDn = {[1 0], [0 1]}
                v = DsDn{1};
                up = geometry.LDdispHS(P(1), P(2), tc.LDxe, tc.LDye, tc.LDa, tc.LDbeta, v(1), v(2), tc.nu);
                um = geometry.LDdispHS(M(1), M(2), tc.LDxe, tc.LDye, tc.LDa, tc.LDbeta, v(1), v(2), tc.nu);
                jump = up - um;
                tc.verifyEqual(norm(jump), 1, 'AbsTol', 1e-5, ...
                    sprintf('|jump| for [Ds Dn] = [%d %d]', v));
                if v(2) == 1
                    tc.verifyEqual(abs(jump*n'), 1, 'AbsTol', 1e-5, 'opening jump is along the normal');
                else
                    tc.verifyEqual(abs(jump*n'), 0, 'AbsTol', 1e-5, 'shear jump is along the element');
                end
            end
        end

        function lineDislocationStressIsHookeOfDisplacement(tc)
            % LDstressHS must equal Hooke's law applied to the gradient of LDdispHS
            G = tc.G; nu = tc.nu; E = 2*G*(1 + nu); lam = 2*G*nu/(1 - 2*nu);
            pts = [20 -25; -8 -3; 4 -14];
            h = 1e-4;
            u = @(x, y) geometry.LDdispHS(x, y, tc.LDxe, tc.LDye, tc.LDa, tc.LDbeta, 1, 0.3, nu);
            for k = 1:size(pts, 1)
                p = pts(k,:);
                dux = (u(p(1)+h, p(2)) - u(p(1)-h, p(2))) / (2*h);
                duy = (u(p(1), p(2)+h) - u(p(1), p(2)-h)) / (2*h);
                exx = dux(1); eyy = duy(2); exy = (duy(1) + dux(2))/2;
                Sfd = [lam*(exx+eyy) + 2*G*exx, lam*(exx+eyy) + 2*G*eyy, 2*G*exy];
                S = geometry.LDstressHS(p(1), p(2), tc.LDxe, tc.LDye, tc.LDa, tc.LDbeta, 1, 0.3, nu, E);
                tc.verifyLessThan(norm(S - Sfd) / norm(S), 1e-6, sprintf('point (%g, %g)', p));
            end
        end

        %% ---------------- triangle strain volume ----------------
        function triangleDisplacementMatchesTanhSinhQuadrature(tc)
            % analytic (Barbot 2018) vs double-exponential numerical integration
            [u2, u3] = computeDisplacementPlaneStrainTriangleShearZone(tc.ObsX2, tc.ObsX3, ...
                tc.TriA, tc.TriB, tc.TriC, 0.3, 1, -0.3, tc.nu);
            [t2, t3] = computeDisplacementPlaneStrainTriangleShearZoneTanhSinh(tc.ObsX2, tc.ObsX3, ...
                tc.TriA, tc.TriB, tc.TriC, 0.3, 1, -0.3, tc.nu);
            scale = max(abs([u2(:); u3(:)]));
            tc.verifyLessThan(max(abs([u2(:)-t2(:); u3(:)-t3(:)])) / scale, 1e-6);
        end

        function triangleDisplacementWithNearlyAlignedEdges(tc)
            % An edge horizontal or vertical up to floating-point noise
            % (a 1e-10 m tilt over 10 km, as PDE Toolbox meshes produce)
            % used to give u2 errors of ~66 per unit strain; the closed
            % form is exact only for exactly aligned edges. Compared with
            % tanh-sinh quadrature around the centroid.
            nu = 0.25; s = (-5:5)'*1e-2;
            cases = {[145e3 156698.7; 140e3 165359; 150e3 165359 + 1e-10], ...   % nearly horizontal base
                     [140e3 160e3; 140e3 + 1e-10 170e3; 150e3 165e3]};             % nearly vertical edge
            for k = 1:numel(cases)
                V = cases{k}; P = mean(V) + s*[1 0];
                for e = {[1 0 -1], [0 1 0]}
                    [u2, u3] = computeDisplacementPlaneStrainTriangleShearZone(P(:,1), P(:,2), V(1,:), V(2,:), V(3,:), e{1}(1), e{1}(2), e{1}(3), nu);
                    [r2, r3] = computeDisplacementPlaneStrainTriangleShearZoneTanhSinh(P(:,1), P(:,2), V(1,:), V(2,:), V(3,:), e{1}(1), e{1}(2), e{1}(3), nu);
                    scale = max(abs([r2; r3]));
                    tc.verifyEqual([u2 u3], [r2 r3], 'AbsTol', 1e-6*scale, sprintf('case %d', k));
                end
            end
        end

        function triangleStressIsSmoothOnEdgeLines(tc)
            % The closed-form displacement has atan(N/D) terms with D = 0 on
            % the lines through the edges of the triangle and of its image
            % above the surface. On those lines (away from the triangle's
            % own edges) the field is smooth. The earlier finite-difference
            % stress (removed) gave stresses ~1e10 too large across the
            % atan jump (structured meshes aligned with the surface put
            % centroids exactly there); kept as a regression case.
            G = 30e3; nu = 0.25;
            A = [0 10e3]; B = [10e3 10e3]; C = [10e3 20e3];     % depth positive
            % on the image edge line A'C' (x3 = -10 km - x2) only, and on
            % the source edge line AC (x3 = 10 km + x2) outside the triangle;
            % with the unit normal of each line
            pts = [-30e3 20e3; 30e3 40e3];
            nrm = [1 1; 1 -1]/sqrt(2);
            for k = 1:size(pts, 1)
                P = pts(k,:); d = 50;   % m: curvature error of the mean ~2e-5
                Q = [P; P + d*nrm(k,:); P - d*nrm(k,:)];
                for e = {[1 0 -1], [0 1 0]}
                    [s22, s23, s33] = computeStressPlaneStrainTriangleShearZone( ...
                        Q(:,1), Q(:,2), A, B, C, e{1}(1), e{1}(2), e{1}(3), G, nu);
                    S = [s22 s23 s33];
                    ref = (S(2,:) + S(3,:))/2;      % smooth field: mean of the two neighbours
                    tc.verifyEqual(S(1,:), ref, 'AbsTol', 1e-4*norm(ref), sprintf('point %d', k));
                end
            end
        end
        %% ---------------- analytic triangle stress ----------------
        function triangleStressAnalyticMatchesDisplacementDerivative(tc)
            % Independent of the generated code: 4th-order central
            % differences of the closed-form displacement (step 1e-3 l,
            % truncation ~1e-12) at points away from the triangle's edges
            % and from the lines where the displacement itself cancels.
            % The last point is inside the triangle, where the eigenstrain
            % is removed: stress = Hooke(grad u - e).
            G = tc.G; nu = tc.nu;
            T = {[0 10e3], [10e3 12e3], [4e3 20e3]};
            X = [-25e3 3e3; 30e3 15e3; 5e3 40e3; -8e3 22e3; 18e3 30e3; 6e3 2e3; 4.5e3 14e3];
            in = inpolygon(X(:,1), X(:,2), [T{1}(1) T{2}(1) T{3}(1)], [T{1}(2) T{2}(2) T{3}(2)]);
            assert(isequal(find(in), size(X,1)), 'only the last point is inside');
            h = 1e-3*norm(T{2} - T{1});
            for e = {[1 0 -1], [0 1 0], [0.3 -0.2 0.5]}
                s = e{1};
                [s22, s23, s33] = computeStressPlaneStrainTriangleShearZone(X(:,1), X(:,2), T{:}, s(1), s(2), s(3), G, nu);
                u = @(dx, dz) displacement(X(:,1) + dx, X(:,2) + dz, T, s, nu);
                d4 = @(f) (-f(2) + 8*f(1) - 8*f(-1) + f(-2))/(12*h);
                g22 = d4(@(k) u(k*h, 0)*[1; 0]);  g32 = d4(@(k) u(k*h, 0)*[0; 1]);
                g23 = d4(@(k) u(0, k*h)*[1; 0]);  g33 = d4(@(k) u(0, k*h)*[0; 1]);
                % elastic strain: displacement gradient minus the eigenstrain inside
                E22 = g22 - in*s(1); E23 = (g23 + g32)/2 - in*s(2); E33 = g33 - in*s(3);
                lambda = 2*G*nu/(1 - 2*nu);
                ref = [lambda*(E22 + E33) + 2*G*E22, 2*G*E23, lambda*(E22 + E33) + 2*G*E33];
                S = [s22 s23 s33];
                tc.verifyLessThan(max(abs(S - ref), [], 2) ./ vecnorm(ref, 2, 2), 1e-7, mat2str(s));
            end
            function U = displacement(x2, x3, T, s, nu)
                [u2, u3] = computeDisplacementPlaneStrainTriangleShearZone(x2, x3, T{:}, s(1), s(2), s(3), nu);
                U = [u2 u3];
            end
        end

        function triangleStressAnalyticSurfaceIsTractionFree(tc)
            % including a triangle with an edge on the surface
            I = eye(3);
            xs = linspace(-60, 80, 141)';
            tris = {{tc.TriA, tc.TriB, tc.TriC}, {[-3 0], [5 0], [1 6]}};
            for t = 1:2
                for k = 1:3
                    [s22, s23, s33] = computeStressPlaneStrainTriangleShearZone(xs, 0*xs, tris{t}{:}, I(k,1), I(k,2), I(k,3), tc.G, tc.nu);
                    ok = abs(xs - 1) > 4.5;       % off the surface edge
                    tc.verifyTrue(all(isfinite([s22(ok); s23(ok); s33(ok)])), sprintf('triangle %d, strain %d: finite', t, k));
                    tc.verifyLessThan(max(abs([s23(ok); s33(ok)])) / max(abs(s22(ok))), 1e-10, ...
                        sprintf('triangle %d, strain %d: s23, s33 at the surface', t, k));
                end
            end
        end

        function triangleStressAnalyticIsReciprocal(tc)
            % Betti with exact cell averages (16 x 16 Gauss product rule):
            % holds to rounding, even for neighbouring triangles
            T1 = {[10 20], [12 20], [11 22]};
            e1 = [0.2 0.7 -0.4]; e2 = [-0.5 0.3 0.1];
            for dx = [2.5 10]
                T2 = {[10 + dx 21], [12 + dx 22], [10.5 + dx 23.5]};
                [P1, w1] = triangleQuadrature(T1); [P2, w2] = triangleQuadrature(T2);
                [s22, s23, s33] = computeStressPlaneStrainTriangleShearZone(P2(:,1), P2(:,2), T1{:}, e1(1), e1(2), e1(3), tc.G, tc.nu);
                w12 = w2'*(s22*e2(1) + 2*s23*e2(2) + s33*e2(3));
                [s22, s23, s33] = computeStressPlaneStrainTriangleShearZone(P1(:,1), P1(:,2), T2{:}, e2(1), e2(2), e2(3), tc.G, tc.nu);
                w21 = w1'*(s22*e1(1) + 2*s23*e1(2) + s33*e1(3));
                tc.verifyEqual(w12, w21, 'RelTol', 1e-9, sprintf('dx = %g', dx));
            end
            function [P, w] = triangleQuadrature(T)
                n = 16; b = (1:n-1)./sqrt(4*(1:n-1).^2 - 1);
                [V, D] = eig(diag(b, 1) + diag(b, -1)); [g, i] = sort(diag(D)); wg = 2*V(1,i)'.^2;
                [U, W] = ndgrid((g + 1)/2, (g + 1)/2); [WU, WW] = ndgrid(wg/2, wg/2);
                l1 = U(:); l2 = W(:).*(1 - U(:));
                A = cell2mat(T'); ar = abs(det([A(2,:) - A(1,:); A(3,:) - A(1,:)]));
                P = A(1,:) + l1*(A(2,:) - A(1,:)) + l2*(A(3,:) - A(1,:));
                w = WU(:).*WW(:).*(1 - U(:))*ar;
            end
        end

        function triangleStressAnalyticIsSmoothOnEdgeLines(tc)
            % On the line through an image edge, and through a horizontal
            % or vertical source edge, the closed form is 0/0 and loses
            % accuracy nearby; the routine must still be continuous there.
            G = tc.G; nu = tc.nu;
            A = [0 10e3]; B = [10e3 10e3]; C = [10e3 20e3];
            % image AC, horizontal AB, vertical BC (below and above), source
            % AC (no special case), image AC far away
            pts = [-30e3 20e3; -20e3 10e3; 10e3 30e3; 10e3 2e3; 30e3 40e3; -300e3 290e3];
            nrm = [1 1; 0 1; 1 0; 1 0; 1 -1; 1 1]; nrm = nrm./vecnorm(nrm, 2, 2);
            for k = 1:size(pts, 1)
                d = [0 1e-6 1e-4]';     % m
                Q = pts(k,:) + d.*nrm(k,:);
                for e = {[1 0 -1], [0 1 0]}
                    [s22, s23, s33] = computeStressPlaneStrainTriangleShearZone(Q(:,1), Q(:,2), A, B, C, e{1}(1), e{1}(2), e{1}(3), G, nu);
                    S = [s22 s23 s33];
                    tc.verifyTrue(all(isfinite(S(:))), sprintf('point %d: finite', k));
                    % the field itself changes by up to ~3e-8 over 1e-4 m
                    tc.verifyEqual(S(2:3,:), repmat(S(1,:), 2, 1), 'AbsTol', 1e-7*norm(S(1,:)), sprintf('point %d', k));
                end
            end
        end

        function triangleStressAcceptsSeveralStrains(tc)
            % [1 x K] strains give [N x K] stresses, equal to K calls
            X = [tc.ObsX2', tc.ObsX3'];
            E = [1 0 0.3; 0 1 -0.2; 0 0 0.5];
            [a, b, c] = computeStressPlaneStrainTriangleShearZone(X(:,1), X(:,2), tc.TriA, tc.TriB, tc.TriC, E(1,:), E(2,:), E(3,:), tc.G, tc.nu);
            tc.verifySize(a, [size(X, 1), 3]);
            for k = 1:3
                [p, q, r] = computeStressPlaneStrainTriangleShearZone(X(:,1), X(:,2), tc.TriA, tc.TriB, tc.TriC, E(1,k), E(2,k), E(3,k), tc.G, tc.nu);
                tc.verifyEqual([a(:,k) b(:,k) c(:,k)], [p q r], 'AbsTol', 1e-12*max(abs([p; q; r])));
            end
        end
    end
end

