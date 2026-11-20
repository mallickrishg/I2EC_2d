classdef (SharedTestFixtures = {I2ECPathFixture}) BoundaryCorrectionTest < matlab.unittest.TestCase
    %BOUNDARYCORRECTIONTEST  The BEM box correction in computeShzStressKernelsBem.
    %
    % computeShzStressKernelsBem adds two parts: the analytic half-space
    % stress of each triangle (geometry.computeShzStressKernels) and the
    % stress of the slip on the box walls that enforces the box conditions
    % (the BEM correction).  Each part on its own is a stress field of the
    % same linear elastic body, so each must satisfy Betti's reciprocal
    % theorem.  The half-space part is checked in KernelsTest; this class
    % checks the correction alone, which the full-kernel check cannot see
    % because the correction is ~1% of the total.
    %
    % Betti weights: for receiver r, response component p (sxx, sxz, szz)
    % and source s, strain component q (exx, exz, ezz), the matrix
    %     B((p,r),(q,s)) = LL(r,s,p,q) * area(r) * w(p),   w = [1 2 1]
    % (the factor 2 counts both shear components in the work) is symmetric.
    %
    % Until 2026-11-19 the correction rows were stacked as (xx, zz, xz)
    % while the half-space rows are (xx, xz, zz), which made the correction
    % ~65% asymmetric on the example mesh.

    properties
        Sm     % small model struct (makeSmallModel)
    end

    methods (TestClassSetup)
        function buildModel(tc)
            folder = tc.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture).Folder;
            tc.Sm = makeSmallModel(folder);
        end
    end

    methods (Test)
        function boxCorrectionIsReciprocal(tc)
            shz = tc.Sm.shz;
            [~, Lhs] = evalc('geometry.computeShzStressKernels(shz, shz)');
            [~, Lbem] = evalc('computeShzStressKernelsBem(shz, shz, tc.Sm.bnd)');
            B = bettiMatrix(Lbem - Lhs, shz.area(:));
            asym = norm(B - B.', 'fro') / norm(B, 'fro');
            tc.verifyLessThan(asym, 1e-2, sprintf('box correction asymmetry %.3g', asym));
        end

        function boxCorrectionIsNotZero(tc)
            % guard: every stress component gets a non-trivial correction,
            % so the reciprocity test above cannot pass on zeros
            shz = tc.Sm.shz;
            [~, Lhs] = evalc('geometry.computeShzStressKernels(shz, shz)');
            [~, Lbem] = evalc('computeShzStressKernelsBem(shz, shz, tc.Sm.bnd)');
            C = Lbem - Lhs;
            for p = 1:3
                tc.verifyGreaterThan(norm(reshape(C(:,:,p,:), [], 1)), ...
                    1e-6 * norm(reshape(Lhs(:,:,p,:), [], 1)), sprintf('component %d', p));
            end
        end
    end
end

function B = bettiMatrix(LL, a)
w = [1 2 1];
N = size(LL, 1);
B = zeros(3*N);
for p = 1:3
    for q = 1:3
        B((p-1)*N+(1:N), (q-1)*N+(1:N)) = LL(:,:,p,q) .* a * w(p);
    end
end
end
