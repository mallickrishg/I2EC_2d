function stresskernel = assembleStressKernel(evl,locked)
% assemble stress kernel for linear ODE system
% from kernel data structure and locked domain

stresskernel = [evl.KK(~locked,~locked),  evl.LK(~locked,:,1),   evl.LK(~locked,:,2);...
                evl.KL(:,~locked,1),      evl.LL(:,:,1,1),       evl.LL(:,:,1,2);...
                evl.KL(:,~locked,2),      evl.LL(:,:,2,1),       evl.LL(:,:,2,2)];

end