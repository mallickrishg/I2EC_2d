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
S = [];
for i=1:3
    for j=1:3
        S.LL{i,j}=zeros(shz.N);
    end
end
% source strain 100;010;001
I=eye(3);
tic
for k=1:shz.N
    for i=1:3
        [s22,s23,s33]=computeStressPlaneStrainTriangleShearZoneFiniteDifference( ...
            shz.xc(:,1),-shz.xc(:,2),...
            shz.A(k,:),shz.B(k,:),shz.C(k,:),...
            I(i,1),I(i,2),I(i,3),mu,nu);
        
        S.LL{1,i}(:,k)=s22(:);
        S.LL{2,i}(:,k)=s23(:);
        S.LL{3,i}(:,k)=s33(:);
    end
end
toc

% construct deviatoric stress kernel and remove positive eigen values

end