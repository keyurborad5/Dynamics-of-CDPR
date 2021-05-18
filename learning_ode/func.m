function [out] = func(t,X)
%FUNC Summary of this function goes here
%   Detailed explanation goes here

k=3;
m=1;
b=2*sqrt(k*m);
k_ps=m/k*10;
k_vs=2*m/k*sqrt(k_ps);
fd=3*rectangularPulse(0,7,t)+10*rectangularPulse(13,15,t);;
F = fd+m*(-k_ps*(X(1)-fd)-k_vs*X(2))/k;
A=[0 1;-k/m -b/(m)]*[X(1);X(2)]+[0;F*k/m];
% x1_dot = X(2)
% x2_dot = (F-k*X(1))/m;
power = F*X(2);
out = [A(1,:)' A(2,:)' power]';
end

