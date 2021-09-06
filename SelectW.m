function [SelectedW,pA] = SelectW(sA,TotalN,Ac)
%
pA = sA;
A1 = sort(pA(pA>0));
A0 = sort(pA(pA<0),'descend');
Iter = length(find(pA~=0))-TotalN;
if nargin == 2
    Ac = 0;
end
for sn = 1:Iter
    if abs(A1(1)+Ac)>abs(A0(1)+Ac)
        pA(pA==A0(1)) = 0;
        Ac = Ac+A0(1);
        if length(A0)>1
            A0(1) = [];
        else
            A0(1) = -10000;
        end
    elseif abs(A1(1)+Ac)<abs(A0(1)+Ac)
        pA(pA==A1(1)) = 0;
        Ac = Ac+A1(1);
        if length(A1)>1
            A1(1) = [];
        else
            A1(1) = 10000;
        end
    end
end
SelectedW = find(pA~=0);
end

