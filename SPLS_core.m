function [B,E,W,U,V,C,D] = SPLS_core(X,Y,H,para)
%PLSN 
%   paraY apply sparse constrain on Y
[n,p]=size(X);
[q]=size(Y,2);
Xh = X;
Yh = Y;
E = ones(n,H); % X scores
W = ones(n,H); % Y scores
U = ones(p,H); % X weights
V = ones(q,H); % Y weights
C = ones(p,H); % X loadings
D = ones(q,H); % Y loadings
convergence = 0.000001;
for h = 1:H
    diff = 10;
    M = Xh'*Yh./norm(Xh'*Yh);
    [Us,~,Vs] = svd(M,0);
    uo = Us(:,1);
    uo = uo./norm(uo);
    vo = Vs(:,1);
    vo = vo./norm(uo);
    while diff > convergence 
        Mvo = M*vo;
        
        %if mean(abs(Mvo))>0.032     for corn
        %    para2 = para*(mean(abs(Mvo))-0.031)*400;
        %else
            para2=para;
        %end
        un = wthresh(Mvo,'s',para2);
        vn = M'*uo; % wthresh(M'*uo,'s',paraY);
        un = un./norm(un);
        vn = vn./norm(vn);
        uo=un;
        vo=vn;          
        diff = max(norm(un - uo),norm(vn - vo)); % Check convergence
    end
    
    eh = Xh*un/(un'*un);
    wh = Yh*vn/(vn'*vn);
    ch = Xh'*eh/(eh'*eh); % weight X
    dh = Yh'*eh/(eh'*eh);
    % rh = Yh'*wh/(wh'*wh); % For canonical mode
    Xh = Xh - eh*ch';
    Yh = Yh - eh*dh';

    U(:,h) = un;
    V(:,h) = vn;       
    E(:,h) = eh;            
    W(:,h) = wh;
    C(:,h) = ch;
    D(:,h) = dh;
end
B= U*((C'*U)^(-1))*D'; % Regression coefficients
end

