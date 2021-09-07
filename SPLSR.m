function [Yf,selectW,B,finalN,OH,RMSECV] = SPLSR(Xc,Xv,Yc,H,para,fold)
%SPLSR 
Ypc = (Yc-mean(Yc,1))./std(Yc);
if nargin > 5
    CVgroup = cvpartition(size(Ypc,1),'KFold',fold); 
    RMSEl = zeros(fold,H);
    for f = 1:fold   
        Vind = CVgroup.test(f);  
        Cind = CVgroup.training(f); 
        RMSE = SPLS_LVs(Xc(Cind,:),Ypc(Cind,:),Xc(Vind,:),Ypc(Vind,:),H,para);  
        RMSEl(f,:) = RMSE;
    end
    RMSECV = mean(RMSEl,1);
    [~,OH] = min(RMSECV);
    B = SPLS_core(Xc,Ypc,OH,para);
    Yf = Xv*B.*std(Yc)+mean(Yc,1);
    %dlmwrite('LVs_SPLS',OH,'-append');
else
    B = SPLS_core(Xc,Ypc,H,para);
    Yf = Xv*B.*std(Yc)+mean(Yc,1);
end
C = size(Ypc,2);
m = size(Xc,2);
selectW = zeros(m,C);
finalN = zeros(1,C);
for c= 1:C
    b = B(:,c);
    n0 = find(b~=0);
    finalN(c) = length(n0);
    selectW(1:finalN(c),c) = n0;
end
%dlmwrite('Ws_SPLS_S14',finalN,'-append');
end

