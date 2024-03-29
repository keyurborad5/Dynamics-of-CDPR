clc
clear all
count=0;
count_x=0;
tmax=0;
tmax_x=0;
fes_points=[];
op_fes_points=[];
tic
%% First of all we will create a outer frame of CDPR 
% Frame of dimension (1mx1mx1m) and its corners are named with suffix A
% supplied coordinates are in world/Base frame
syms A;
A = [0,0,0;
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
syms B;
% B=0.01*[-5,5,-5;        %4
%         5,-5,5;         %6
%         -5,-5,-5;       %1
%         5,5,5;          %7
%          5,-5,-5;       %2
%          -5,5,5;        %8
%          5,5,-5;        %3
%          -5,-5,5];      %5
     %%%%%%%%%%%%%%%%%%%%%%%%%As given in papers
     B=0.01*[-10,5,-5;        %4
        10,-5,5;         %6
        -10,-5,-5;       %1
        10,5,5;          %7
         10,-5,-5;       %2
         -10,5,5;        %8
         10,5,-5;        %3
         -10,-5,5];      %5
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% B = 0.01*[-5,-5,-5;
%     5,-5,-5;
%     5,5,-5;
%     -5,5,-5;
%     -5,-5,5;
%     5,-5,5;
%     5,5,5;
%     -5,5,5];
    
%% Now position of COM and Rotation matrix wrt base frame
syms R P phi psi theta;
Px=0:0.025:1;
% Px=0.5;
% Py=0.51;
% Pz=0.5;
Py=0:0.025:1;
Pz=0:0.1:1;
loop_no=0;
 
[X,Y,Z]=meshgrid(Px,Py,Pz);
for i_1=1:length(Px)
    for j_1=1:length(Py)
        for k_1=1:length(Pz)
            P=[X(i_1,j_1,k_1);Y(i_1,j_1,k_1);Z(i_1,j_1,k_1)];
            loop_no=loop_no+1
% P = [0.3;0.3;0.4]; %position of COM of object in base frame
%fixed angle rotation x(psi),y(theta),z(phi)
%RXYZ = Rz(?)*Ry(?)*Rx(?)
psi=0;phi=0*pi/180;theta=0;
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
    L(:,i)=a(i,:)'-(P+R*B(i,:)');
    Lm(i)=norm(a(i,:)'-(P+R*B(i,:)'));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% L=[L ([P(1);P(2);0]-P)]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
L';
Lm;
%% Now to plot the cable from object end to frame end
O = zeros(3,8);
for i=1:length(O)
    O(:,i) = (P+R*B(i,:)');
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
cross_prod(:,i) = cross(R*B(i,:)',unit_v(i,:)');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% cross_prod = [cross_prod cross(P,unit_v(9,:)')]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
W = [unit_v';cross_prod];
%r_k = rank(W);
%% Wrench applied to the platform due to gravity Wg
m_p = 1;%mass of platform as 1KG
MS_p = R*m_p*[0;0;0];
g=[0;0;-9.81];
Wg = [m_p*eye(3);
        [0 -MS_p(3) MS_p(2);
        MS_p(3) 0 -MS_p(1);
        -MS_p(2) MS_p(1) 0]]*g;
%% Finding tension in the string
% W*t+We+Wg=0
% here t=pinv(W)*(-Wg)Moore-Penrose pseudoinverse
% t = pinv(W)*(-Wg-[0;0;9;4.5;-4.5;0]); %+ null(W,'r')*[2;1]
% t = pinv(W)*[0,0,9.81*4,0,0,0]';
% norm(t);
% 
% if any(t<4) || any(t>80)
%     count=count+1;
% else
%     if max(t)>tmax
%         tmax=max(t);
%     end
%     fes_points=[fes_points ; P']; %saving feasible points in an array
%     figure(2)
%     plot3(P(1),P(2),P(3),".r")
%     axis([-0.1 1.5 -0.1 1.5 -0.2 1.1])
%     xlabel('x-axis');
%     ylabel('y-axis');
%     zlabel('z-axis');
%     grid on
%     hold on
% end
%         end
%     end
% end
%% Quadratic programming optimisation for the cable force
fmin=4;
fmax=80;
f_ref=(fmax*2)/2;
H=eye(8);
f=-f_ref*ones(1,8);
lb=fmin*ones(8,1);
ub=fmax*ones(8,1);
Aeq=W;
Beq=[0,0,9.81*4,0,0,0]';
options = optimset('Display', 'off');
[x fval exitflag] = quadprog(H,f,[],[],Aeq,Beq,lb,ub,[],options);

if exitflag ~=1
    count_x=count_x+1;
else
    if max(x)>tmax_x
        tmax_x=max(x);
    end
       op_fes_points=[op_fes_points ; P'];
%     figure(3)
%     plot3(P(1),P(2),P(3),".r")
%     axis([-0.1 1.5 -0.1 1.5 -0.2 1.1])
%     xlabel('x-axis');
%     ylabel('y-axis');
%     zlabel('z-axis');
%     grid on
%     hold on
end
        end
    end
end
toc

%********************************
%% For plotting saved feasible points (should run this by coping in command window, donot uncomment here)
% for i=1:length(fes_points)
%figure(1)
% plot3(fes_points(i,1),fes_points(i,2),fes_points(i,3),".r")
%      axis([-0.1 1.5 -0.1 1.5 -0.2 1.1])
%      xlabel('x-axis');
%      ylabel('y-axis');
%      zlabel('z-axis');
%      grid on
%      hold on
% end
%*****************************