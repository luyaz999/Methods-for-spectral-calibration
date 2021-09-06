function [Hpar] = GetParameters()
%GETPARAMEERS 
Hpar.Iter = 80; %
Hpar.ro = 1.1; % 
Hpar.mu = 0.0002; %
Hpar.lm = 0.01;
Hpar.lm2 = 0.03;
Hpar.Ncomp = 11;
Hpar.lm = 0.001; %0.002 for SCWSR
Hpar.lm2 = 0.0001;
%Hpar.Ncomp = 27; %35;
%Hpar.NcompS = 28;
Hpar.paraS = 0.07;
Hpar.funlm = @(x) (1);
Hpar.m = 3;
end

