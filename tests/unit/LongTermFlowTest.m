classdef (SharedTestFixtures = {I2ECPathFixture}) LongTermFlowTest < matlab.unittest.TestCase
    %LONGTERMFLOWTEST  Corner-flow loading: computeCornerflowVelocityfield,
    % compute_cornerflow_velocityfield (old), getStrainratesLongterm,
    % computeLongTermVelocity.
    %
    % Corner flow is checked against its own boundary conditions
    % (arc: v = 0 at the surface, v = slab velocity on the slab; ocean:
    % v = (1,0) on the incoming plate, slab velocity on the slab) and
    % incompressibility.  Velocities are normalised by the plate rate.
    %
    % Momentum balance is checked independently of how the coefficients
    % were found (cornerFlowStokesResidual): for a corner flow, power-law
    % Stokes reduces to h'' + k^2 h = 0 with h = F^(1/n), k^2 = (2n-1)/n^2.
    % The n = 3 checks run over a range of dips, including shallow ones
    % where the coefficients are large and nearly cancel.

    properties (Constant)
        Dip = 20*pi/180
        R = [0.5; 1; 2]
    end

    properties (TestParameter)
        dipDeg = {3, 5, 8, 14, 20, 30, 45, 60, 80}
        n = {1, 3}
    end

    properties
        Ex
    end

    methods (TestClassSetup)
        function loadExample(tc)
            tc.Ex = loadExampleModel(false);
        end
    end

    methods (Test)
        function newtonianCornerFlowBoundaryConditions(tc)
            tc.checkBoundaryConditions(1, 1e-10);
        end

        function powerLawCornerFlowBoundaryConditions(tc, dipDeg)
            tc.checkBoundaryConditions(3, 1e-8, dipDeg*pi/180);
        end

        function cornerFlowSolvesStokes(tc, dipDeg, n)
            % both regions: momentum balance for exponent n, velocity
            % independent of r
            d = dipDeg*pi/180;
            regions = {linspace(-d, 0, 401), linspace(-pi + 0.01, -d, 401)};
            names = {'wedge', 'ocean'};
            for col = 1:2
                f = @(x, z) tc.cornerFlowColumn(d, x, z, n, col);
                [stokes, radial] = cornerFlowStokesResidual(f, regions{col}(2:end-1), n);
                tc.verifyLessThan(stokes, 1e-3, sprintf('%s: momentum balance', names{col}));
                tc.verifyLessThan(radial, 1e-8, sprintf('%s: velocity depends on r', names{col}));
            end
        end

        function stokesCheckRejectsWrongRheology(tc)
            % guard for the check above: the Newtonian wedge flow is not a
            % power-law (n = 3) solution, and vice versa
            d = tc.Dip; th = linspace(-d, 0, 401); th = th(2:end-1);
            f1 = @(x, z) tc.cornerFlowColumn(d, x, z, 1, 1);
            f3 = @(x, z) tc.cornerFlowColumn(d, x, z, 3, 1);
            tc.verifyGreaterThan(cornerFlowStokesResidual(f1, th, 3), 0.05);
            tc.verifyGreaterThan(cornerFlowStokesResidual(f3, th, 1), 0.05);
        end

        function powerLawCornerFlowIsQuiet(tc)
            clear computeCornerflowVelocityfield   % empty the coefficient cache
            out = evalc('computeCornerflowVelocityfield(tc.Dip, 1, -0.1, 3);');
            tc.verifyEmpty(strtrim(out), 'console output');
        end

        function powerLawCoefficientsAreSolvedOncePerDip(tc)
            % getStrainratesLongterm evaluates one triangle per call
            clear computeCornerflowVelocityfield
            computeCornerflowVelocityfield(tc.Dip, 1, -0.1, 3);
            t = tic;
            for k = 1:200
                computeCornerflowVelocityfield(tc.Dip, rand, -rand, 3);
            end
            tc.verifyLessThan(toc(t), 1, '200 calls at the same dip');
        end

        function cornerFlowIsIncompressible(tc)
            [X, Z] = tc.wedgeGrid();
            arc = atan2(Z, X) > -tc.Dip;
            h = 1e-6;
            for n = [1 3]
                f = @(x, z) tc.quietCornerFlow(x, z, n);
                [vxp, ~] = f(X+h, Z); [vxm, ~] = f(X-h, Z);
                [~, vzp] = f(X, Z+h); [~, vzm] = f(X, Z-h);
                div = (vxp - vxm)/(2*h) + (vzp - vzm)/(2*h);
                tc.verifyLessThan(max(abs(div(arc,1))), 1e-4, sprintf('arc, n = %d', n));
                tc.verifyLessThan(max(abs(div(~arc,2))), 1e-4, sprintf('ocean, n = %d', n));
            end
        end

        function oldAndNewNewtonianCornerFlowAgree(tc)
            % compute_cornerflow_velocityfield (used by computeLongTermVelocity)
            % vs computeCornerflowVelocityfield(n = 1) (used by getStrainratesLongterm)
            [X, Z] = tc.wedgeGrid();
            [vx1, vz1] = computeCornerflowVelocityfield(tc.Dip, X, Z, 1);
            [vx0, vz0] = compute_cornerflow_velocityfield(tc.Dip, X, Z);
            tc.verifyEqual(vx0, vx1, 'AbsTol', 1e-10*max(abs(vx1(:))));
            tc.verifyEqual(vz0, vz1, 'AbsTol', 1e-10*max(abs(vz1(:))));
        end

        function invalidExponentIsRejected(tc)
            tc.verifyError(@() computeCornerflowVelocityfield(tc.Dip, 1, -1, 2), ?MException);
        end

        function longTermStrainRatesAreFiniteAndDeviatoric(tc)
            shz = tc.Ex.shz;
            [e22, e23] = getStrainratesLongterm(shz, tc.Ex.rcv.dip(1)*pi/180, [0 20e3], [-140e3 35e3], 1);
            tc.verifySize(e22, [shz.N 1]);
            tc.verifySize(e23, [shz.N 1]);
            tc.verifyTrue(all(isfinite([e22; e23])));
            tc.verifyGreaterThan(max(abs([e22; e23])), 0);
        end

        function longTermStrainRatesMatchVelocityGradient(tc, n)
            % the planar fit over the triangle vertices must reproduce the
            % finite-difference gradient of the same corner-flow field at
            % the centroid (same translations as getStrainratesLongterm)
            shz = tc.Ex.shz; dip = tc.Ex.rcv.dip(1)*pi/180;
            oc = [0 20e3 - 1]; co = [-140e3 35e3 - 1];
            [e22, e23] = getStrainratesLongterm(shz, dip, [0 20e3], [-140e3 35e3], n);
            cont = shz.xc(:,2) + shz.xc(:,1)*tan(dip) >= 0;
            h = 1;   % m
            E22 = zeros(shz.N,1); E23 = E22;
            for i = 1:shz.N
                if cont(i), t = co; col = 1; else, t = oc; col = 2; end
                x = shz.xc(i,1) + t(1); z = shz.xc(i,2) + t(2);
                [a, b] = computeCornerflowVelocityfield(dip, [x+h; x-h; x; x], [z; z; z+h; z-h], n);
                ux = (a(1,col) - a(2,col))/(2*h); uz = (a(3,col) - a(4,col))/(2*h);
                vx = (b(1,col) - b(2,col))/(2*h); vz = (b(3,col) - b(4,col))/(2*h);
                E22(i) = 0.5*(ux - vz); E23(i) = 0.5*(uz + vx);
            end
            relerr = hypot(e22 - E22, e23 - E23) ./ hypot(E22, E23);
            tc.verifyLessThan(median(relerr), 0.05, 'median relative error');
        end

        function longTermVelocityUsesSameCornerFlowAsStrainRates(tc)
            % computeLongTermVelocity (old corner flow, polygon regions) vs the
            % corner flow behind getStrainratesLongterm, in the oceanic mantle
            rcv = tc.Ex.rcv; shz = tc.Ex.shz; bnd = tc.Ex.bnd;
            dipdeg = rcv.dip(1);
            X = [150e3; 250e3; 300e3; 100e3];
            Z = [-150e3; -200e3; -300e3; -250e3];   % below the slab
            [Vx, Vz] = computeLongTermVelocity(X, Z, 1, 20e3, dipdeg, shz, bnd, [0 20e3], [-140e3 35e3]);
            [a, b] = computeCornerflowVelocityfield(dipdeg*pi/180, X + 0, Z + 20e3 - 1, 1);
            tc.verifyEqual(Vx, a(:,2), 'AbsTol', 1e-10);
            tc.verifyEqual(Vz, b(:,2), 'AbsTol', 1e-10);
        end

        function strainRatesWithoutExponentDefaultToNewtonian(tc)
            % n is optional and defaults to 1 (older scripts call with 4 arguments)
            shz = tc.Ex.shz; dip = tc.Ex.rcv.dip(1)*pi/180;
            [e22, e23] = getStrainratesLongterm(shz, dip, [0 20e3], [-140e3 35e3]);
            [r22, r23] = getStrainratesLongterm(shz, dip, [0 20e3], [-140e3 35e3], 1);
            tc.verifyEqual([e22 e23], [r22 r23]);
        end
    end

    methods (Access = private)
        function checkBoundaryConditions(tc, n, tol, d)
            if nargin < 4, d = tc.Dip; end
            r = tc.R; U = [cos(d) -sin(d)];
            [vx, vz] = computeCornerflowVelocityfield(d, r, 0*r, n);                         % arc surface
            tc.verifyLessThan(max(abs([vx(:,1); vz(:,1)])), tol, 'arc: v = 0 at surface');
            [vx, vz] = computeCornerflowVelocityfield(d, r*cos(d), -r*sin(d), n);            % arc slab
            tc.verifyLessThan(max(abs([vx(:,1)-U(1); vz(:,1)-U(2)])), tol, 'arc: slab velocity');
            [vx, vz] = computeCornerflowVelocityfield(d, -r, 0*r - 1e-12, n);                % incoming plate
            tc.verifyLessThan(max(abs([vx(:,2)-1; vz(:,2)])), tol, 'ocean: plate velocity');
            [vx, vz] = computeCornerflowVelocityfield(d, r*cos(d), -r*sin(d) - 1e-12, n);    % ocean slab
            tc.verifyLessThan(max(abs([vx(:,2)-U(1); vz(:,2)-U(2)])), tol, 'ocean: slab velocity');
        end

        function [vx, vz] = quietCornerFlow(tc, x, z, n)
            [vx, vz] = computeCornerflowVelocityfield(tc.Dip, x, z, n);
        end

        function [vx, vz] = cornerFlowColumn(~, d, x, z, n, col)
            [Vx, Vz] = computeCornerflowVelocityfield(d, x, z, n);
            vx = Vx(:, col); vz = Vz(:, col);
        end

        function [X, Z] = wedgeGrid(~)
            [X, Z] = meshgrid(linspace(0.1, 2, 15), -linspace(0.05, 1.5, 15));
            X = X(:); Z = Z(:);
        end
    end
end
