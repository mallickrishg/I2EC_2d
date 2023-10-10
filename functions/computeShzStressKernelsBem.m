function LL = computeShzStressKernelsBem(src,shz,boundary)

% initialize stress kernels
LL1 = zeros(shz.N,src.N);
LL2 = zeros(shz.N,src.N);
LL3 = zeros(shz.N,src.N);

% this are 3x3 stress kernels
LL = zeros(shz.N,shz.N,3,3);

% source strain 100;010;001
I = eye(3);

tic
disp('Computing shear zone - shear zone stress kernels')
% convert all depths to positive numbers
xc = shz.xc; xc(:,2) = -xc(:,2);
A = shz.A;A(:,2) = -A(:,2);
B = shz.B;B(:,2) = -B(:,2);
C = shz.C;C(:,2) = -C(:,2);

for i = 1:3
    % each iteration of 'i' goes through each eigen strain source
    % i = 1 corresponds to e22 source
    % i = 2 corresponds to e23 source
    % i = 3 corresponds to e33 source
    parfor k = 1:src.N
        [s22,s23,s33] = computeStressPlaneStrainTriangleShearZoneFiniteDifference( ...
            xc(:,1),xc(:,2),...
            A(k,:),B(k,:),C(k,:),...
            I(i,1),I(i,2),I(i,3),...
            mu,nu);
        
        % need to work with 2-d matrices because MATLAB doesn't like 3-d or
        % 4-d matrices inside parfor
        LL1(:,k) = s22(:);
        LL2(:,k) = s23(:);
        LL3(:,k) = s33(:);

    end
    LL(:,:,1,i) = LL1;
    LL(:,:,2,i) = LL2;
    LL(:,:,3,i) = LL3;
end

end