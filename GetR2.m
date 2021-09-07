function [R2] = GetR2(Yf,Y)
%GETR2 
C = size(Y,2);
R2 = zeros(1,C);
for ic = 1:C
R2(ic) = 1 - norm(Y(:,ic)-Yf(:,ic))^2/norm(Y(:,ic) - mean(Y(:,ic)))^2;
end
end

