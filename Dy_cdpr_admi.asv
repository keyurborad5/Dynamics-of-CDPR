clc
clear all
tic
no=0;
yes=0;

params.B=0.01*[-10,15,-5;        %4
        10,-15,5;         %6
        -10,-15,-5;       %1
        10,15,5;          %7
         10,-15,-5;       %2
         -10,15,5;        %8
         10,15,-5;        %3
         -10,-15,5];      %5
%*****************************
params.a=[0,0.05,1;
    0.05,0,1;
    0.95,0,1;
    1,0.05,1;
    1,0.95,1;
    0.95,1,1;
    0.05,1,1;
    0,0.95,1];
%***************************
%% Wrench applied to the platform due to gravity Wg
mp = 100;%mass of platform as 1KG
params.mp=mp;
R=eye(3,3);
MSp = R*mp*[0;0;0];
params.MSp = MSp;
g=[0;0;-9.81];
params.Wg = [mp*eye(3);
        [0 -MSp(3) MSp(2);
        MSp(3) 0 -MSp(1);
        -MSp(2) MSp(1) 0]]*g;
box_l=0.2;
box_b=0.4;
box_h=0.1;
I_xx=mp*(box_h^2+box_b^2)/12;
params.I_xx=I_xx;
I_yy=mp*(box_l^2+box_h^2)/12;
params.I_yy=I_yy;
I_zz=mp*(box_l^2+box_b^2)/12;
params.I_zz=I_zz;
params.Ig = [I_xx 0 0;
        0 I_yy 0;
        0 0 I_zz];


%MSp = mp*R*[0;0;0];
params.MSp_hat = [0 -MSp(3) MSp(2);
            MSp(3) 0 -MSp(1);
            -MSp(2) MSp(1) 0];
% Ip=R*Ig*R'- MSp_hat^2/mp;
% IIp=[mp*eye(3,3) -MSp_hat;
%     MSp_hat Ip];

%% for virtual object
mp_v=1;
l_v=0.2;
b_v=0.4;
h_v=0.1;
I_xx_v=mp_v*(h_v^2+b_v^2)/12;
I_yy_v=mp_v*(l_v^2+h_v^2)/12;
I_zz_v=mp_v*(l_v^2+b_v^2)/12;
Ig_v = [I_xx_v 0 0;
        0 I_yy_v 0;
        0 0 I_zz_v];

M_v=[mp_v*eye(3) zeros(3);
    zeros(3) Ig_v];
params.inv_M_v=inv(M_v);
%*************************************
t_span = 0:0.1:50;
inti_x = [0.5;0.5;0.75;0;0;0;0;0;0;0;0;0];%0.5;0.5;0.5;0;0;0;0;0;0;0;0;0];
[t X]= ode45(@(t,X) Dycdpradmi(t,X,params), t_span, inti_x);

%% Optimising the string forces
y=0.5*rectangularPulse(3,7,t)+1*rectangularPulse(13,18,t)-0.5*rectangularPulse(25,35,t)-1*rectangularPulse(43,46,t);
%y=0.5*rectangularPulse(0,7,t)+1*rectangularPulse(13,18,t);
fes_x=sparse(length(t),6);
syms R P phi psi theta;
% Px=0:0.1:1.5;
% % Px=0.5;
% % Py=0.51;
% % Pz=0.5;
% Py=0:0.1:1.5;
% Pz=0:0.1:1;
% loop_no=0;
%  
% [X,Y,Z]=meshgrid(Px,Py,Pz);
% for i_1=1:length(Px)
%     for j_1=1:length(Py)
%         for k_1=1:length(Pz)
%             P=[X(i_1,j_1,k_1);Y(i_1,j_1,k_1);Z(i_1,j_1,k_1)];
%             loop_no=loop_no+1
for iter=1:length(t)
P = [X(iter,1);X(iter,2);X(iter,3)]; %position of COM of object in base frame
%fixed angle rotation x(psi),y(theta),z(phi)
%RXYZ = Rz(?)*Ry(?)*Rx(?)
psi=X(iter,4);phi=X(iter,6);theta=X(iter,5);
R = [cos(phi)*cos(theta) -sin(phi)*cos(psi)+cos(phi)*sin(theta)*sin(psi) sin(phi)*sin(psi)+cos(phi)*sin(theta)*cos(psi);
    sin(phi)*cos(theta) cos(phi)*cos(psi)+sin(phi)*sin(theta)*sin(psi) -cos(phi)*sin(psi)+sin(phi)*sin(theta)*cos(psi);
    -sin(theta) cos(theta)*sin(psi) cos(theta)*cos(psi)]; %Rotation of frame of Object wrt Base frame
%% Now vector representing the direction of the string and length of string is givrn as below
L = zeros(3,8);
Lm = zeros(8,1);
% a=[A(5,:);
%     A(5,:);
%     A(6,:);
%     A(6,:);
%     A(7,:);
%     A(7,:);
%     A(8,:);
%     A(8,:)];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
a=[0,0.05,1;
    0.05,0,1;
    0.95,0,1;
    1,0.05,1;
    1,0.95,1;
    0.95,1,1;
    0.05,1,1;
    0,0.95,1];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:length(L)
    L(:,i)=a(i,:)'-(P+R*params.B(i,:)');
    Lm(i)=norm(a(i,:)'-(P+R*params.B(i,:)'));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% L=[L ([P(1);P(2);0]-P)]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
L';
Lm;
%% Now to plot the cable from object end to frame end
O = zeros(3,8);
for i=1:length(O)
    O(:,i) = (P+R*params.B(i,:)');
end
% plot3([O(1,3) O(1,5)],[O(2,3) O(2,5)],[O(3,3) O(3,5)], ...
%     [O(1,5) O(1,7)],[O(2,5) O(2,7)],[O(3,5) O(3,7)], ...
%     [O(1,7) O(1,1)],[O(2,7) O(2,1)],[O(3,7) O(3,1)], ...
%     [O(1,1) O(1,3)],[O(2,1) O(2,3)],[O(3,1) O(3,3)], ...
%     [O(1,3) O(1,8)],[O(2,3) O(2,8)],[O(3,3) O(3,8)], ...
%     [O(1,8) O(1,2)],[O(2,8) O(2,2)],[O(3,8) O(3,2)], ...
%     [O(1,2) O(1,4)],[O(2,2) O(2,4)],[O(3,2) O(3,4)], ...
%     [O(1,4) O(1,6)],[O(2,4) O(2,6)],[O(3,4) O(3,6)], ...
%     [O(1,6) O(1,8)],[O(2,6) O(2,8)],[O(3,6) O(3,8)], ...
%     [O(1,6) O(1,1)],[O(2,6) O(2,1)],[O(3,6) O(3,1)], ...
%     [O(1,4) O(1,7)],[O(2,4) O(2,7)],[O(3,4) O(3,7)], ...
%     [O(1,2) O(1,5)],[O(2,2) O(2,5)],[O(3,2) O(3,5)], ...
%     'LineWidth',2)
%     hold on
% plot3(P(1),P(2),P(3),'o','MarkerFaceColor','red','MarkerSize',10)
% axis on
% xlabel('x')
% ylabel('y')
% zlabel('z')
% hold on
% for i=1:8
%     plot3([a(i,1) O(1,i)],[a(i,2) O(2,i)],[a(i,3) O(3,i)],'b','LineWidth',1)
%     hold on
%     grid on
%     axis([0 1 -0 1 -0 1])
% end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % plot3([P(1) P(1)],[P(2) P(2)],[P(3) 0],'g','LineWidth',1)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% hold off 
%% Static and kinematic model
% Wrench matrics
unit_v=zeros(8,3);
for i=1:8
unit_v(i,:)=L(:,i)'/norm(L(:,i));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% unit_v = [unit_v; L(:,9)'/norm(L(:,9))];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:8
cross_prod(:,i) = cross(R*params.B(i,:)',unit_v(i,:)');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% cross_prod = [cross_prod cross(P,unit_v(9,:)')]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
W = [unit_v';cross_prod];

%% Quadratic programming optimisation for the cable force
fmin=4;
fmax=80;
f_ref=(fmax+fmin)/2;
H=eye(8);
f=-f_ref*ones(1,8);
lb=fmin*ones(8,1);
ub=fmax*ones(8,1);
Aeq=W;
Beq=[-y(i),0,9.81*3,0,0,0]';
options = optimset('Display', 'off');
[x fval exitflag] = quadprog(H,f,[],[],Aeq,Beq,lb,ub,[],options);

if exitflag ~= 1
    fes_x(iter,:)=fes_x(iter-1,:);
    no=no+1;
else
    fes_x(iter,:)=X(iter,1:6);
    yes=yes+1;
end
end



figure(1)
%plot(t,X(:,1)/k,'r',t,X(:,2)/k,'--g',t,X(:,1),'b',t,y,'k')
plot(t,X(:,1),'r',t,X(:,2),'g',t,X(:,3),'--c',t,y,'k')
grid on
title('Virtual_position')

figure(2)
%plot(t,X(:,1)/k,'r',t,X(:,2)/k,'--g',t,X(:,1),'b',t,y,'k')
plot(t,fes_x(:,1),'r',t,fes_x(:,2),'g',t,fes_x(:,3),'--c',t,y,'k')
grid on
title('Feasible_position')
% figure(2)
% %plot(t,X(:,1)/k,'r',t,X(:,2)/k,'--g',t,X(:,1),'b',t,y,'k')
% plot(t,X(:,13),'m',t,X(:,14),'g',t,X(:,15),'--y',t,y,'k')
% grid on
% title('actcual_position')

%% Now position of COM and Rotation matrix wrt base frame
%syms phi psi theta;
% P = [X(:,13)';X(:,14)';X(:,15)']; %position of COM of object in base frame
% %fixed angle rotation x(psi),y(theta),z(phi)
% %RXYZ = Rz(phi)*Ry(theta)*Rx(psi)
% psi=X(:,16)*pi/180;phi=X(:,18)*pi/180;theta=X(:,17)*pi/180; %provide angles in radians
% for i=1:length(t)
% R(:,:,i) = [cos(phi(i))*cos(theta(i)) -sin(phi(i))*cos(psi(i))+cos(phi(i))*sin(theta(i))*sin(psi(i)) sin(phi(i))*sin(psi(i))+cos(phi(i))*sin(theta(i))*cos(psi(i));
%     sin(phi(i))*cos(theta(i)) cos(phi(i))*cos(psi(i))+sin(phi(i))*sin(theta(i))*sin(psi(i)) -cos(phi(i))*sin(psi(i))+sin(phi(i))*sin(theta(i))*cos(psi(i));
%     -sin(theta(i)) cos(theta(i))*sin(psi(i)) cos(theta(i))*cos(psi(i))]; %Rotation of frame of Object wrt Base frame
% end
% R=real(R);
% %% Now to plot the cable from object end to frame end
% O = zeros(3,8,length(P));
% for j=1:length(P)
%     for i=1:8
%         O(:,i,j) = (P(:,j)+R(:,:,j)*params.B(i,:)');
%     end
%     
%     figure(3)
%     plot3([O(1,3,j) O(1,5,j)],[O(2,3,j) O(2,5,j)],[O(3,3,j) O(3,5,j)], ...
%     [O(1,5,j) O(1,7,j)],[O(2,5,j) O(2,7,j)],[O(3,5,j) O(3,7,j)], ...
%     [O(1,7,j) O(1,1,j)],[O(2,7,j) O(2,1,j)],[O(3,7,j) O(3,1,j)], ...
%     [O(1,1,j) O(1,3,j)],[O(2,1,j) O(2,3,j)],[O(3,1,j) O(3,3,j)], ...
%     [O(1,3,j) O(1,8,j)],[O(2,3,j) O(2,8,j)],[O(3,3,j) O(3,8,j)], ...
%     [O(1,8,j) O(1,2,j)],[O(2,8,j) O(2,2,j)],[O(3,8,j) O(3,2,j)], ...
%     [O(1,2,j) O(1,4,j)],[O(2,2,j) O(2,4,j)],[O(3,2,j) O(3,4,j)], ...
%     [O(1,4,j) O(1,6,j)],[O(2,4,j) O(2,6,j)],[O(3,4,j) O(3,6,j)], ...
%     [O(1,6,j) O(1,8,j)],[O(2,6,j) O(2,8,j)],[O(3,6,j) O(3,8,j)], ...
%     [O(1,6,j) O(1,1,j)],[O(2,6,j) O(2,1,j)],[O(3,6,j) O(3,1,j)], ...
%     [O(1,4,j) O(1,7,j)],[O(2,4,j) O(2,7,j)],[O(3,4,j) O(3,7,j)], ...
%     [O(1,2,j) O(1,5,j)],[O(2,2,j) O(2,5,j)],[O(3,2,j) O(3,5,j)], ...
%     'LineWidth',2)
%     hold on
%     plot3(P(1,j),P(2,j),P(3,j),'o','MarkerFaceColor','red','MarkerSize',10)
%     axis on
%     xlabel('x')
%     ylabel('y')
%     zlabel('z')
%     hold on
%     for i=1:8
%         plot3([params.a(i,1) O(1,i,j)],[params.a(i,2) O(2,i,j)],[params.a(i,3) O(3,i,j)],'b','LineWidth',1)
%         hold on
%         grid on
%         axis([0 1 0 1 -0.3 1.5])
%     end
%     hold off
% end
toc

%% FUNCTION 
function [out] = Dycdpradmi(t,X,params)
% mp_v=10;
% l_v=0.2;
% b_v=0.4;
% h_v=0.1;
% I_xx_v=mp_v*(h_v^2+b_v^2)/12;
% I_yy_v=mp_v*(l_v^2+h_v^2)/12;
% I_zz_v=mp_v*(l_v^2+b_v^2)/12;
% Ig_v = [I_xx_v 0 0;
%         0 I_yy_v 0;
%         0 0 I_zz_v];
% 
% M_v=[mp_v*eye(3) zeros(3);
%     zeros(3) Ig_v];
% inv_M_v=inv(M_v);
inv_M_v=params.inv_M_v;
pos_v=[X(1) X(2) X(3)];
norm_pos_v = norm(pos_v-[0.5 0.5 0.75]);
if norm_pos_v < 0.1
    K_v=5*eye(6);
else
    K_v=0*eye(6);
end
B_v =10*eye(6);

fmea=[0.5*rectangularPulse(0,7,t)+1*rectangularPulse(13,18,t)-0.5*rectangularPulse(25,35,t)-1*rectangularPulse(43,46,t);0;0;0;0;0];
%fmea=[0.5*rectangularPulse(0,7,t)+1*rectangularPulse(13,18,t);0;0;0;0;0];
fdes=[0;0;0;0;0;0]; %it will have -mg in z direction here i have not consider it in fmea or fdes

States_v = [X(1);X(2);X(3);X(4);X(5);X(6);X(7);X(8);X(9);X(10);X(11);X(12)];
out_v=[zeros(6) eye(6);
        -inv_M_v*K_v -inv_M_v*B_v]*States_v+[0;0;0;0;0;0;inv_M_v*[fmea-fdes+K_v*[0.5;0.5;0.75;0;0;0]]];

%position control
% kps=30;
% bps=80;
% F=kps*[X(1)-X(13);X(2)-X(14);X(3)-X(15);X(4)-X(16);X(5)-X(17);X(6)-X(18);]+...
%     bps*[X(7)-X(19);X(8)-X(20);X(9)-X(21);X(10)-X(22);X(11)-X(23);X(12)-X(24)];
% %actual system dynamics
% %fixed angle rotation x(psi),y(theta),z(phi)
% psi=X(16)*pi/180;phi=X(18)*pi/180;theta=X(17)*pi/180; %provide angles in radians
% R = [cos(phi)*cos(theta) -sin(phi)*cos(psi)+cos(phi)*sin(theta)*sin(psi) sin(phi)*sin(psi)+cos(phi)*sin(theta)*cos(psi);
%     sin(phi)*cos(theta) cos(phi)*cos(psi)+sin(phi)*sin(theta)*sin(psi) -cos(phi)*sin(psi)+sin(phi)*sin(theta)*cos(psi);
%     -sin(theta) cos(theta)*sin(psi) cos(theta)*cos(psi)]; %Rotation of frame of Object wrt Base frame
% R=real(R);
% mp=params.mp;
% Ip=R*params.Ig*R'- params.MSp_hat^2/mp;
% IIp=[mp*eye(3,3) -params.MSp_hat;
%     params.MSp_hat Ip];
% invIIp=eye(6)/IIp;
% B = [[params.MSp(2) params.MSp(3) 0; params.MSp(1) 0 params.MSp(3);0 params.MSp(1) params.MSp(2)];
%     [0 0 params.I_zz-params.I_yy; 0 params.I_xx-params.I_zz 0; params.I_yy-params.I_xx 0 0]];
% C = [[0 -params.MSp(1) -params.MSp(1); -params.MSp(2) 0 -params.MSp(2);-params.MSp(3) -params.MSp(3) 0];
%      zeros(3)];
% K = 0;
% d = [0;0;-mp*9.81;0;0;0];
% States_a = [X(19);X(20);X(21);X(22);X(23);X(24);X(22)*X(23);X(22)*X(24);X(23)*X(24);X(22)^2;X(23)^2;X(24)^2];
% out_a=[eye(6) zeros(6);
%         zeros(6) -invIIp*B -invIIp*C]*States_a+[0;0;0;0;0;0;invIIp*[F]];
out = [out_v']';
end