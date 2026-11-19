classdef (SharedTestFixtures = {I2ECPathFixture}) GreensFunctionsTest < matlab.unittest.TestCase
    %GREENSFUNCTIONSTEST  Elementary half-space solutions (lowest level of I2EC).
    %
    % Line dislocations (geometry.LDdispHS / geometry.LDstressHS, from
    % CutAndDisplace) and triangular strain volumes (Barbot 2018:
    % computeDisplacementPlaneStrainTriangleShearZone*, computeStress*).
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

        function triangleStressFiniteDifferenceSurfaceIsTractionFree(tc)
            % the stress routine used by every shear-zone kernel
            I = eye(3);
            for k = 1:3
                [s22, s23, s33] = computeStressPlaneStrainTriangleShearZoneFiniteDifference( ...
                    [-30 0 25 60], [0 0 0 0], tc.TriA, tc.TriB, tc.TriC, I(k,1), I(k,2), I(k,3), tc.G, tc.nu);
                scale = max(abs(s22));
                tc.verifyLessThan(max(abs([s23 s33])) / scale, 1e-5, ...
                    sprintf('s23, s33 at the surface, source strain %d', k));
            end
        end

        function triangleStressFiniteDifferenceIsReciprocal(tc)
            % Betti: work of eigenstrain e1 (triangle 1) on the stress of e2
            % (triangle 2) equals the reverse.  Small triangles, far apart,
            % so centroid collocation is accurate.
            T1 = {[10 20], [12 20], [11 22]};
            T2 = {[40 30], [42 31], [40.5 32.5]};
            a1 = polyarea([T1{1}(1) T1{2}(1) T1{3}(1)], [T1{1}(2) T1{2}(2) T1{3}(2)]);
            a2 = polyarea([T2{1}(1) T2{2}(1) T2{3}(1)], [T2{1}(2) T2{2}(2) T2{3}(2)]);
            c1 = mean(cell2mat(T1'));
            c2 = mean(cell2mat(T2'));
            e1 = [0.2 0.7 -0.4];  % [e22 e23 e33]
            e2 = [-0.5 0.3 0.1];
            [s22, s23, s33] = computeStressPlaneStrainTriangleShearZoneFiniteDifference( ...
                c2(1), c2(2), T1{:}, e1(1), e1(2), e1(3), tc.G, tc.nu);
            w12 = a2 * (s22*e2(1) + 2*s23*e2(2) + s33*e2(3));
            [s22, s23, s33] = computeStressPlaneStrainTriangleShearZoneFiniteDifference( ...
                c1(1), c1(2), T2{:}, e2(1), e2(2), e2(3), tc.G, tc.nu);
            w21 = a1 * (s22*e1(1) + 2*s23*e1(2) + s33*e1(3));
            tc.verifyEqual(w12, w21, 'RelTol', 1e-3);
        end
    end

    methods (Test, TestTags = {'KnownIssue'})
        function triangleStressAnalyticMatchesFiniteDifference(tc)
            % computeStressPlaneStrainTriangleShearZone.m (2,238 lines, not
            % called anywhere) should agree with the finite-difference
            % version that the kernels use.  Evaluated point by point
            % because the analytic routine rejects vector input (next test).
            I = eye(3);
            maxRel = 0;
            for k = 1:3
                for p = 1:numel(tc.ObsX2)
                    [a1, a2, a3] = computeStressPlaneStrainTriangleShearZone(tc.ObsX2(p), tc.ObsX3(p), ...
                        tc.TriA, tc.TriB, tc.TriC, I(k,1), I(k,2), I(k,3), tc.G, tc.nu);
                    [f1, f2, f3] = computeStressPlaneStrainTriangleShearZoneFiniteDifference(tc.ObsX2(p), tc.ObsX3(p), ...
                        tc.TriA, tc.TriB, tc.TriC, I(k,1), I(k,2), I(k,3), tc.G, tc.nu);
                    maxRel = max(maxRel, norm([a1 a2 a3] - [f1 f2 f3]) / norm([f1 f2 f3]));
                end
            end
            tc.verifyLessThan(maxRel, 1e-3, ...
                'analytic and finite-difference triangle stresses disagree (see u2d2 overwrite at line 84)');
        end

        function triangleStressAnalyticAcceptsVectorInput(tc)
            I = eye(3);
            f = @() computeStressPlaneStrainTriangleShearZone(tc.ObsX2(:), tc.ObsX3(:), ...
                tc.TriA, tc.TriB, tc.TriC, I(1,1), I(1,2), I(1,3), tc.G, tc.nu);
            tc.verifyWarningFree(f, 'Omega(x2,x3)*e22 uses matrix product, fails for vectors');
        end
    end
end
