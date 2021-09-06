function [SelectedW,pA,MaxN,A,LossFun] = SWS(X,Y,Hpar,TotalN)
%Sparse coefficients wavelength selection and regression method selects feature wavelengths of X for predicting Y.
%  The inputs are the dictionary (D), the LR HSI image (X), the number of atoms in dictionary (K),
%    number of iteration (J) and hyper parameters (lm,mu,ro).
%  The output are the coefficient matrix (B) and recorded fun of B (fun).
%I/O: [SelectedW,A,funB,TotalN] = NSWS(X,Y,Hpar);
Y = (Y-mean(Y,1))./std(Y); % Std for soil
[~, L] = size(X);
Iter = Hpar.Iter;
lm = Hpar.lm;
mu = Hpar.mu;
ro = Hpar.ro;
C = size(Y, 2);
A = zeros(L,C)+0.15; % Sparse matrix
U = A; % Lagrangian multiplier
LossFun = zeros(1,Iter);
XTX = X' * X;
XTY = X' * Y;
I = eye(L); % Identity matrix
MaxN = zeros(C,Iter);
%SelectedW = zeros(TotalN,C);
pA = zeros(L,C);
for  j  =  1:Iter
    S = (XTX + 2 * mu * I)^(-1)*(XTY + 2 * mu * (A - U/(2 * mu)));
    A = wthresh(S + U /(2 * mu),'s',lm /(2 * mu));
    U  = U + mu*(S - A);
    mu = mu * ro;
    %LossFun(j) = 0.5*sum(sum((Y-X*A).^2)) + lm*sum(sum(abs(A)));
    [nr,nc] = find(A~=0);
    for ic = 1:C
        MaxN(ic,j) = length(nr(nc==ic));
    end
end
if nargin > 3
    SelectedW = zeros(max(TotalN),C);
    sA = A.*repmat(mean(X,1)',1,C); % Wavelength scores matrix
    Yf = sum(sA,1);
    Ac = mean(Y,1)-Yf;
    for ic = 1:C
        [SelectedW(1:min(sum(sA(:,ic)~=0),TotalN(ic)),ic),pA(:,ic)] = SelectW(sA(:,ic),TotalN(ic),Ac(ic));     
        %[~,IndexS] = sort(abs(A(:,ic)),'descend'); 
        %SelectedW(1:TotalN(ic),ic) = IndexS(1:TotalN(ic));
    end
else
    SelectedW = zeros(max(MaxN(:,Iter)),C);
    for ic = 1:C
        SelectedW(1:MaxN(ic,Iter),ic) = nr(nc==ic);
    end
end
end

