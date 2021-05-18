clc
clear all
tic

%% First of all we will create a outer frame of CDPR 
% Frame of dimension (1mx1mx1m) and its corners are named with suffix A
% supplied coordinates are in world/Base frame
params.A = [0,0,0;
    1,0,0;
    1,1,0;
    0,1,0;
    0,0,1;
    1,0,1;
    1,1,1;
    0,1,1];
%% Secondly will create an Object (cube) which we want to manipulate of dimension (0.1m,0.1m,0.1m)
%supplied coordinate will be in objects frame
%object frame is located at the COM of the object
% B=0.01*[-5,5,-5;        %4
%         5,-5,5;         %6
%         -5,-5,-5;       %1
%         5,5,5;          %7
%          5,-5,-5;       %2
%          -5,5,5;        %8
%          5,5,-5;        %3
%          -5,-5,5];      %5
     %%%%%%%%%%%%%%%%%%%%%%%%%As given in papers
     params.B=0.01*[-10,15,-5;        %4
        10,-15,5;         %6
        -10,-15,-5;       %1
        10,15,5;          %7
         10,-15,-5;       %2
         -10,15,5;        %8
         10,15,-5;        %3
         -10,-15,5];      %5
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% B = 0.01*[-10,-20,-5;
%     10,-20,-5;
%     10,20,-5;
%     -10,20,-5;
%     -10,-20,5;
%     10,-20,5;
%     10,20,5;
%     -10,20,5];
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
mp = 10;%mass of platform as 1KG
params.mp=mp;
R=eye(3,3);
MSp = R*mp*[0;0;0];
g=[0;0;-9.81];
params.Wg = [mp*eye(3);
        [0 -MSp(3) MSp(2);
        MSp(3) 0 -MSp(1);
        -MSp(2) MSp(1) 0]]*g;
box_l=0.2;
box_b=0.4;
box_h=0.1;
I_xx=mp*(box_h^2+box_b^2)/12;
I_yy=mp*(box_l^2+box_h^2)/12;
I_zz=mp*(box_l^2+box_b^2)/12;
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

%****************************
t_span = 0:0.1:30;
inti_x = [0.5;0.5;0.5;0;0;0;0;0;0;0;0;0];
% m=1
% k=3
[t X]= ode45(@(t,X) dynaCDPR(t,X,params), t_span, inti_x);
% for i = 1:length(t)
%     Total_energy(i) = 0.5*m/k*X(i,2)^2 + 0.5*k/k*X(i,1)^2+0*2*sqrt(2)*X(i,1)*X(i,2);
% end
% Total_energy;
% X(:,3);
% diff=Total_energy'-X(:,3);
y=3*rectangularPulse(3,7,t)-10*rectangularPulse(13,15,t);
figure(1)
%plot(t,X(:,1)/k,'r',t,X(:,2)/k,'--g',t,X(:,1),'b',t,y,'k')
plot(t,X(:,1),'r',t,X(:,2),'--b',t,X(:,3),'g',t,y,'k')
grid on
title('position')
%% Now position of COM and Rotation matrix wrt base frame
%syms phi psi theta;
P = [X(:,1)';X(:,2)';X(:,3)']; %position of COM of object in base frame
%fixed angle rotation x(psi),y(theta),z(phi)
%RXYZ = Rz(phi)*Ry(theta)*Rx(psi)
psi=X(:,4)*pi/180;phi=X(:,6)*pi/180;theta=X(:,5)*pi/180; %provide angles in radians
for i=1:length(t)
R(:,:,i) = [cos(phi(i))*cos(theta(i)) -sin(phi(i))*cos(psi(i))+cos(phi(i))*sin(theta(i))*sin(psi(i)) sin(phi(i))*sin(psi(i))+cos(phi(i))*sin(theta(i))*cos(psi(i));
    sin(phi(i))*cos(theta(i)) cos(phi(i))*cos(psi(i))+sin(phi(i))*sin(theta(i))*sin(psi(i)) -cos(phi(i))*sin(psi(i))+sin(phi(i))*sin(theta(i))*cos(psi(i));
    -sin(theta(i)) cos(theta(i))*sin(psi(i)) cos(theta(i))*cos(psi(i))]; %Rotation of frame of Object wrt Base frame
end
R=real(R);
%% Now to plot the cable from object end to frame end
O = zeros(3,8,length(P));
for j=1:length(P)
    for i=1:8
        O(:,i,j) = (P(:,j)+R(:,:,j)*params.B(i,:)');
    end
    
    figure(2)
    plot3([O(1,3,j) O(1,5,j)],[O(2,3,j) O(2,5,j)],[O(3,3,j) O(3,5,j)], ...
    [O(1,5,j) O(1,7,j)],[O(2,5,j) O(2,7,j)],[O(3,5,j) O(3,7,j)], ...
    [O(1,7,j) O(1,1,j)],[O(2,7,j) O(2,1,j)],[O(3,7,j) O(3,1,j)], ...
    [O(1,1,j) O(1,3,j)],[O(2,1,j) O(2,3,j)],[O(3,1,j) O(3,3,j)], ...
    [O(1,3,j) O(1,8,j)],[O(2,3,j) O(2,8,j)],[O(3,3,j) O(3,8,j)], ...
    [O(1,8,j) O(1,2,j)],[O(2,8,j) O(2,2,j)],[O(3,8,j) O(3,2,j)], ...
    [O(1,2,j) O(1,4,j)],[O(2,2,j) O(2,4,j)],[O(3,2,j) O(3,4,j)], ...
    [O(1,4,j) O(1,6,j)],[O(2,4,j) O(2,6,j)],[O(3,4,j) O(3,6,j)], ...
    [O(1,6,j) O(1,8,j)],[O(2,6,j) O(2,8,j)],[O(3,6,j) O(3,8,j)], ...
    [O(1,6,j) O(1,1,j)],[O(2,6,j) O(2,1,j)],[O(3,6,j) O(3,1,j)], ...
    [O(1,4,j) O(1,7,j)],[O(2,4,j) O(2,7,j)],[O(3,4,j) O(3,7,j)], ...
    [O(1,2,j) O(1,5,j)],[O(2,2,j) O(2,5,j)],[O(3,2,j) O(3,5,j)], ...
    'LineWidth',2)
    hold on
    plot3(P(1,j),P(2,j),P(3,j),'o','MarkerFaceColor','red','MarkerSize',10)
    axis on
    xlabel('x')
    ylabel('y')
    zlabel('z')
    hold on
    for i=1:8
        plot3([params.a(i,1) O(1,i,j)],[params.a(i,2) O(2,i,j)],[params.a(i,3) O(3,i,j)],'b','LineWidth',1)
        hold on
        grid on
        axis([0 1 0 1 -0.3 1.5])
    end
    hold off
end

% figure(2)
% plot(t,X(:,2))
% title('velocity')
% figure(3)
% plot(t,Total_energy- X(:,3)')
% title('error')
toc

function [out] = dynaCDPR(t,X,params)
%FUNC Summary of this function goes here
%   Detailed explanation goes here
%% First of all we will create a outer frame of CDPR 
% Frame of dimension (1mx1mx1m) and its corners are named with suffix A
% supplied coordinates are in world/Base frame
%t
% A = [0,0,0;
%     1,0,0;
%     1,1,0;
%     0,1,0;
%     0,0,1;
%     1,0,1;
%     1,1,1;
%     0,1,1];
%% Secondly will create an Object (cube) which we want to manipulate of dimension (0.1m,0.1m,0.1m)
%supplied coordinate will be in objects frame
%object frame is located at the COM of the object
% B=0.01*[-5,5,-5;        %4
%         5,-5,5;         %6
%         -5,-5,-5;       %1
%         5,5,5;          %7
%          5,-5,-5;       %2
%          -5,5,5;        %8
%          5,5,-5;        %3
%          -5,-5,5];      %5
     %%%%%%%%%%%%%%%%%%%%%%%%%As given in papers
%      B=0.01*[-10,20,-5;        %4
%         10,-20,5;         %6
%         -10,-20,-5;       %1
%         10,20,5;          %7
%          10,-20,-5;       %2
%          -10,20,5;        %8
%          10,20,-5;        %3
%          -10,-20,5];      %5
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% B = 0.01*[-10,-20,-5;
%     10,-20,-5;
%     10,20,-5;
%     -10,20,-5;
%     -10,-20,5;
%     10,-20,5;
%     10,20,5;
%     -10,20,5];
    
%% Now position of COM and Rotation matrix wrt base frame
%syms phi psi theta;
P = [X(1);X(2);X(3)]; %position of COM of object in base frame
%fixed angle rotation x(psi),y(theta),z(phi)
%RXYZ = Rz(phi)*Ry(theta)*Rx(psi)
psi=X(4)*pi/180;phi=X(6)*pi/180;theta=X(5)*pi/180; %provide angles in radians
R = [cos(phi)*cos(theta) -sin(phi)*cos(psi)+cos(phi)*sin(theta)*sin(psi) sin(phi)*sin(psi)+cos(phi)*sin(theta)*cos(psi);
    sin(phi)*cos(theta) cos(phi)*cos(psi)+sin(phi)*sin(theta)*sin(psi) -cos(phi)*sin(psi)+sin(phi)*sin(theta)*cos(psi);
    -sin(theta) cos(theta)*sin(psi) cos(theta)*cos(psi)]; %Rotation of frame of Object wrt Base frame
R=real(R);
%% Now vector representing the direction of the string and length of string is givrn as below
L = zeros(3,8);
%Lm = zeros(8,1);
% a=[A(1,:);
%    A(5,:);
%    A(2,:);
%    A(6,:);
%     A(3,:);
%     A(7,:);
%     A(4,:);
%     A(8,:)];
%*****************************
% a=[0,0.05,1;
%     0.05,0,1;
%     0.95,0,1;
%     1,0.05,1;
%     1,0.95,1;
%     0.95,1,1;
%     0.05,1,1;
%     0,0.95,1];
%***************************
%**************************
% a=[A(4,:);
%     A(6,:);
%     A(1,:);
%     A(7,:);
%     A(2,:);
%     A(8,:);
%     A(3,:);
%     A(5,:)];

%************************
for i=1:length(L)
    L(:,i)=params.a(i,:)'-(P+R*params.B(i,:)');
    %Lm(i)=norm(params.a(i,:)'-(P+R*params.B(i,:)'));
end
%L';
%Lm;
% %% Now to plot the cable from object end to frame end
% O = zeros(3,8);
% for i=1:length(O)
%     O(:,i) = (P+R*params.B(i,:)');
% end
% figure(1)
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
%     plot3([params.a(i,1) O(1,i)],[params.a(i,2) O(2,i)],[params.a(i,3) O(3,i)],'b','LineWidth',1)
%     hold on
%     grid on
%     axis([-0.5 1.5 -0.5 1.5 -0.3 1.5])
% end
% hold off 
%% Static and kinematic model
% Wrench matrics
unit_v=zeros(8,3);
for i=1:8
unit_v(i,:)=L(:,i)'/norm(L(:,i));
end
for i=1:8
% cross_prod(:,i) = cross(O(:,i),unit_v(i,:)');
cross_prod(:,i) = cross(R*params.B(i,:)',unit_v(i,:)');
end
W = [unit_v';cross_prod];
r_k = rank(W);
%% Wrench applied to the platform due to gravity Wg
mp = params.mp;%mass of platform as 1KG
% MSp = R*mp*[0;0;0];
% g=[0;0;-9.81];
% Wg = [mp*eye(3);
%         [0 -MSp(3) MSp(2);
%         MSp(3) 0 -MSp(1);
%         -MSp(2) MSp(1) 0]]*g;
%% Finding tension in the string
% W*t+We+Wg=0
% here t=pinv(W)*(-Wg)Moore-Penrose pseudoinverse
fd=3*rectangularPulse(3,7,t)-10*rectangularPulse(13,15,t);
%tsn = pinv(W)*(-params.Wg-[-fd;0;0;0;0;0]) %+ null(W,'r')*[2;1]
%norm(t);
% figure(3)
% plot(t,tsn(1),'-r',t,tsn(2),'-b',t,tsn(3),'-g',t,tsn(4),'-k',t,tsn(5),'-c',t,tsn(6),'-m',t,tsn(7),'-y',t,tsn(8),'--r')
% grid on
% title('tensions')
% hold on
%mp=2;
% box_l=0.2;
% box_b=0.4;
% box_h=0.1;
% I_xx=mp*(box_h^2+box_b^2)/12;
% I_yy=mp*(box_l^2+box_h^2)/12;
% I_zz=mp*(box_l^2+box_b^2)/12;
% Ig = [I_xx 0 0;
%         0 I_yy 0;
%         0 0 I_zz];
% 
% %R=eye(3,3);
% %MSp = mp*R*[0;0;0];
% MSp_hat = [0 -MSp(3) MSp(2);
%             MSp(3) 0 -MSp(1);
%             -MSp(2) MSp(1) 0];
Ip=R*params.Ig*R'- params.MSp_hat^2/mp;
IIp=[mp*eye(3,3) -params.MSp_hat;
    params.MSp_hat Ip];
invIIp=eye(6)/IIp;


%k=0*eye(6);
%m=1;
b=20*eye(6);
%k_ps=m/k*10;
%k_vs=2*m/k*sqrt(k_ps);
%fd=3*rectangularPulse(0,7,t)+10*rectangularPulse(13,15,t);;
%F = fd+m*(-k_ps*(X(1)-fd)-k_vs*X(2))/k;
Fini=[0;0;-mp*9.81;0;0;0];
%Fmea = W*tsn;
Fmea = (-params.Wg-[-fd;0;0;0;0;0]);
OUTPUT=[zeros(6) eye(6);zeros(6) -invIIp*b]*[X(1);X(2);X(3);X(4);X(5);X(6);
                                    X(7);X(8);X(9);X(10);X(11);X(12)]+[0;0;0;0;0;0;invIIp*(Fmea+Fini)];
% x1_dot = X(2)
% x2_dot = (F-k*X(1))/m;
%power = F*X(2);
out = [OUTPUT];
end



