clc
clear all
tic
t_span = 0:0.1:50;
params.c=0;
params.ini_x_v=1;
params.x_0=params.ini_x_v;
inti_x = [params.ini_x_v;0;1;0;0];
m=10;
k=3; % actually k=0 here just mentioned here for calc 'b' 
[t X]= ode45(@(t,X) one_d_admi_fun(t,X,params), t_span, inti_x);
% for i = 1:length(t)
%     Total_energy(i) = 0.5*m/k*X(i,2)^2 + 0.5*k/k*X(i,1)^2+0*2*sqrt(2)*X(i,1)*X(i,2);
% end
% Total_energy;
% X(:,3);
% diff=Total_energy'-X(:,3);
%y=3*rectangularPulse(0,7,t)+10*rectangularPulse(13,15,t)+5*rectangularPulse(25,30,t)-2*rectangularPulse(40,45,t);
y=-3*rectangularPulse(0,7,t)-10*rectangularPulse(13,15,t)-2*rectangularPulse(30,35,t);
figure(1)
%plot(t,X(:,1)/k,'r',t,X(:,2)/k,'--g',t,X(:,1),'b',t,y,'k')
plot(t,X(:,1),'b',t,X(:,3),'r',t,y,'k')
grid on
title('position')
% figure(2)
% plot(t,X(:,2))
% title('velocity')
% figure(3)
% plot(t,Total_energy- X(:,3)')
% title('error')
toc
function [out] = one_d_admi_fun(t,X,params)
%FUNC Summary of this function goes here
%   Detailed explanation goes here
%virtual object model
% syms c ini_x_v
% x_0=0;
% if ini_x_v==X(1)
%     c=c+1;
% else
%     ini_x_v=X(1);
%     c=0;
% end
% c
% X(1)
% if c>=10
%     x_0=X(1);
% else
% end

if X(1) < 2+1 && X(1) > -2+1
    k_v=1.5; %actually k=0 but here written to calc 'b'
else
    k_v=0;
end
m_v=5;
b_v=8;

%position control
kps=30;
bps=80;

%actual system dynamics
M = 50;
B = 0;
K = 0;

f=-3*rectangularPulse(0,7,t)-10*rectangularPulse(13,15,t)-2*rectangularPulse(30,35,t);
fd=0;


tau=kps*(X(1)-X(3))+bps*(X(2)-X(4));
MAT=[0 1 0 0 ;-k_v/m_v -b_v/m_v 0 0;0 0 0 1;0 0 -K/M -B/M];
STATES=[X(1);X(2);X(3);X(4)];
EXT=[0;(f-fd+1*k_v)/m_v;0 ;(tau-f)/M];
A=MAT*STATES+EXT;
% x1_dot = X(2)
% x2_dot = (F-k*X(1))/m;
power = tau*X(2);
out = [A(1,:)' A(2,:)' A(3,:)' A(4,:)' power]';
end