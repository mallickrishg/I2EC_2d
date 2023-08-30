function evl = compute_all_stresskernels(rcv,shz)

mu = rcv.earthModel.G;
nu = rcv.earthModel.nu;

% create data structure to hold all kernels
evl = [];

%% compute fault-fault interaction kernels
[K,~] = computetractionkernels(rcv,rcv);
evl.rcv.KK = K;

%% compute shz-shz interaction kernels
% initialize stress kernels
LL1 = zeros(shz.N,shz.N);
LL2 = zeros(shz.N,shz.N);
LL3 = zeros(shz.N,shz.N);

% source strain 100;010;001
I=eye(3);
tic

xc = shz.xc;
A = shz.A;
B = shz.B;
C = shz.C;

for i=1:3
    for k=1:shz.N
        [s22,s23,s33] = computeStressPlaneStrainTriangleShearZoneFiniteDifference( ...
            xc(:,1),-xc(:,2),...
            A(k,:),B(k,:),C(k,:),...
            I(i,1),I(i,2),I(i,3),...
            mu,nu);

        LL1(:,k) = s22(:);
        LL2(:,k) = s23(:);
        LL3(:,k) = s33(:);

        % LL(:,k,1,i)=s22(:);
        % LL(:,k,2,i)=s23(:);
        % LL(:,k,3,i)=s33(:);

    end
    evl.LL(:,:,1,i) = LL1;
    evl.LL(:,:,2,i) = LL2;
    evl.LL(:,:,3,i) = LL3;
end

toc

% construct deviatoric stress kernel and remove positive eigen values

end