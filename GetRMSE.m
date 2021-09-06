function [RMSE] = GetRMSE(Yf,Y)
%GETRMSE 
C = size(Y,2);
RMSE = zeros(1,C);
for ic = 1:C
    RMSE(ic) = sqrt(mean(mean((Y(:,ic)-Yf(:,ic)).^2)));
end
end

