clc
clear all
tic
t_span = 0:2:20;
inti_x = [0;0;0];
m=1
k=3
[t X]= ode45(@(t,X) func(t,X), t_span, inti_x);
% for i = 1:length(t)
%     Total_energy(i) = 0.5*m/k*X(i,2)^2 + 0.5*k/k*X(i,1)^2+0*2*sqrt(2)*X(i,1)*X(i,2);
% end
% Total_energy;
% X(:,3);
% diff=Total_energy'-X(:,3);
y=3*rectangularPulse(0,7,t)+10*rectangularPulse(13,15,t);
figure(1)
%plot(t,X(:,1)/k,'r',t,X(:,2)/k,'--g',t,X(:,1),'b',t,y,'k')
plot(t,X(:,1)/k,'r',t,y,'k')
grid on
title('position')
% figure(2)
% plot(t,X(:,2))
% title('velocity')
% figure(3)
% plot(t,Total_energy- X(:,3)')
% title('error')
toc