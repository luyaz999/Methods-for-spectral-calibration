function [RMSE] = SPLS_LVs(Xc,Yc,Xv,Yv,H,para)
%SPLS_LVS 
[B,~,~,U,~,C,D] = SPLS_core(Xc,Yc,H,para);
n = size(Yc,2);
RMSE = zeros(H,1);
for h = 1:H-1
    Bh = U(:,1:h)*((C(:,1:h)'*U(:,1:h))^(-1))*D(:,1:h)';
    RMSE(h,:) = sqrt(mean(mean((Yv-Xv*Bh).^2))); %GetRMSE(Xv*Bh,Yv);
end
RMSE(H,:) = sqrt(mean(mean((Yv-Xv*Bh).^2))); %GetRMSE(Xv*B,Yv);
end

