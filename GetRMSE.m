function [RMSE] = GetRMSE(Yf,Y)
%GETRMSE 此处显示有关此函数的摘要
%   此处显示详细说明
C = size(Y,2);
RMSE = zeros(1,C);
for ic = 1:C
    RMSE(ic) = sqrt(mean(mean((Y(:,ic)-Yf(:,ic)).^2)));
end
end

