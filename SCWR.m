function [fYv,SelectedWr,finalN,Ar] = SCWR(Xv,Xc,Yc,Hpar,TotalN)
%SWSR 
%   TotalN: Number of wavelengths

C = size(Yc, 2);
fYv = zeros(size(Xv,1),C);
mYc = mean(Yc,1);
sYc = std(Yc,1); %std for soil % ones(1,C);for others
if nargin > 4
    %Hpar.lm = Hpar.lm*5;
    Ar = zeros(size(Xv,2),C);
    AddN = repmat(Hpar.m,5,C); % 5 for potato
    sN = zeros(1,C); 
    [SelectedWr,~,~,A] = SWS(Xc,Yc,Hpar,TotalN+AddN);  
    finalN = sum(A~=0,1);      
    for ic = 1:C 
        if finalN(ic) <= TotalN(ic)     % Full wavelength regression            
            fYv(:,ic) = Xv*A(:,ic).*sYc(ic)+mYc(ic); 
            Ar(:,ic) = A(:,ic);
        else         
            Hpar.lm = Hpar.lm2*Hpar.funlm(TotalN(ic));
            %Hpar.ro = 1.05; % For corn
            while AddN(ic) >= sN(ic) && AddN(ic) >= 0                   
                AddN(ic) = AddN(ic) - 1;
                iSW = SelectedWr(SelectedWr(:,ic)~=0,ic);    
                [iSWr,~,~,iAr,~] = SWS(Xc(:,iSW),Yc(:,ic),Hpar,TotalN(ic)+AddN(ic));                               
                iSW1 = iSWr(iSWr~=0);    
                SelectedWr(1:length(iSW1),ic) = iSW(iSW1);   
                SelectedWr(length(iSW1)+1:end,ic) = 0;                             
                sN(ic) = sum(iAr==0);                    
            end          
            fYv(:,ic) = Xv(:,iSW)*iAr.*sYc(ic)+mYc(ic);
            finalN(ic) = sum(iAr~=0,1); 
            Ar(iSW,ic) = iAr;
        end   
    end   
else
    %Hpar.lm = Hpar.lm2;
    [SelectedWr,~,~,Ar,~] = SWS(Xc,Yc,Hpar); % Full wavelength regression    
    fYv = Xv*Ar.*sYc+mYc;
    finalN = sum(Ar~=0,1);
     %Record automatically selected wavelengths
    %dlmwrite('SfinalN01',finalN,'-append');
end
end

