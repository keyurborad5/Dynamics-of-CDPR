clc
clear all
tic
yes=0;
no=0;
%% Outer Frame
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
a=[0,0.05,1;
    0.05,0,1;
    0.95,0,1;
    1,0.05,1;
    1,0.95,1;
    0.95,1,1;
    0.05,1,1;
    0,0.95,1];
params.B=0.01*[-10,15,-5;        %4
        10,-15,5;         %6
        -10,-15,-5;       %1
        10,15,5;          %7
         10,-15,-5;       %2
         -10,15,5;        %8
         10,15,-5;        %3
         -10,-15,5];      %5

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
inv_M_v=inv(M_v);
B_v =10*eye(6);
%% Constants of Quadratic programming
fmin=4;
fmax=80;
f_ref=(fmax+fmin)/2;
H=eye(8);
f_op=-f_ref*ones(1,8);
lb=fmin*ones(8,1);
ub=fmax*ones(8,1);
%% defining applied force
time=0:0.1:60;
y=[0.5*rectangularPulse(0,7,time)+1*rectangularPulse(13,18,time)-0.5*rectangularPulse(25,35,time)-1*rectangularPulse(43,48,time)-0.4*rectangularPulse(51,55,time)];
fdes=[0;0;0;0;0;0]; %it will have -mg in z direction here i have not consider it in fmea or fdes
%% define function handle
%X=[disp;velo] ==> disp=X(1,:),velo=X(2,:)
%i made function handle with so many variables because i need to change its
%value while in a loop of solving RK45
f=@(t,X,K_v,x_ref,fmea)[...
    [zeros(6) eye(6);
        -inv_M_v*K_v -inv_M_v*B_v]*[X(1);X(2);X(3);X(4);X(5);X(6);X(7);X(8);X(9);X(10);X(11);X(12)]+[0;0;0;0;0;0;inv_M_v*[fmea-fdes+K_v*[x_ref]]]];%x_ref is column[x ;y; z; a; b; t]
%% initial conditions
t(1)=0;
X(:,1)=[0.5;0.5;0.5;0;0;0;0;0;0;0;0;0];
%% step size
h=0.1;
tfinal=60;
N=ceil(tfinal/h);
x_ref=X(1:6,1);
count=0;
%% RK45
for i=1:N
    %This conditional statements deals with being stiff and compliant
    norm_pos_v = norm(X(1:3,i)-x_ref(1:3));
    if norm_pos_v < 0.1
       K_v=5*eye(6);%behave stiff
       
    else
        K_v=0*eye(6);%behave compliant
        if norm(X(1:3,i)-X(1:3,i-1))<1e-3 %checks if the position is constant for some time
            count=count+1;
            if count >=3/h % that some time is decided by this no. of counts
                x_ref=X(1:6,i); %changes the reference around which it needed to be stiff
            end
        else
            count=0;
        end
    end
    %update time
    t(i+1)=t(i)+h;
    fmea=[0.5*rectangularPulse(0,7,t(i))+1*rectangularPulse(13,18,t(i))-0.5*rectangularPulse(25,35,t(i))-1*rectangularPulse(43,48,t(i))-0.4*rectangularPulse(51,55,t(i));0;0;0;0;0];
    
    %upsade of X
    k1=f(t(i)    ,X(:,i)       ,K_v,x_ref,fmea);
    k2=f(t(i)+h/2,X(:,i)+k1*h/2,K_v,x_ref,fmea);
    k3=f(t(i)+h/2,X(:,i)+k2*h/2,K_v,x_ref,fmea);
    k4=f(t(i)+h  ,X(:,i)+k3*h  ,K_v,x_ref,fmea);
    X(:,i+1)=X(:,i)+(h/6)*(k1+2*k2+2*k3+k4);
    
    %% for workspace feasibility
    P = [X(1,i+1);X(2,i+1);X(3,i+1)]; %position of COM of object in base frame
    %fixed angle rotation x(psi),y(theta),z(phi)
    %RXYZ = Rz(?)*Ry(?)*Rx(?)
    psi=X(4,i+1);phi=X(6,i+1);theta=X(5,i+1);
    R = [cos(phi)*cos(theta) -sin(phi)*cos(psi)+cos(phi)*sin(theta)*sin(psi) sin(phi)*sin(psi)+cos(phi)*sin(theta)*cos(psi);
        sin(phi)*cos(theta) cos(phi)*cos(psi)+sin(phi)*sin(theta)*sin(psi) -cos(phi)*sin(psi)+sin(phi)*sin(theta)*cos(psi);
        -sin(theta) cos(theta)*sin(psi) cos(theta)*cos(psi)]; %Rotation of frame of Object wrt Base frame
    
    L = zeros(3,8);
    Lm = zeros(8,1);
    for j=1:length(L)
        L(:,j)=a(j,:)'-(P+R*params.B(j,:)');
        Lm(j)=norm(a(j,:)'-(P+R*params.B(j,:)'));
    end
    
    %% Static and kinematic model
    % Wrench matrics
    unit_v=zeros(8,3);
    for j=1:8
        unit_v(j,:)=L(:,j)'/norm(L(:,j));
%     end
%     for j=1:8
        cross_prod(:,j) = cross(R*params.B(j,:)',unit_v(j,:)');
    end
    W = [unit_v';cross_prod];
    %% Quadratic programming optimisation for the cable force
    
    Aeq=W;
    %-y means in + dir of x ie (>0.5): 0.6,0.7,0.8..
    %y means in - dir of x ie (<0.5):0.4,0.3,0.2..
    Beq=[0,0,9.81*4,0,0,0]';
    options = optimset('Display', 'off');
    [cab_f fval exitflag] = quadprog(H,f_op,[],[],Aeq,Beq,lb,ub,[],options);
    
    if exitflag ~= 1
        X(:,i+1)=X(:,i-1);
        no=no+1;
    else
        %X(:,i+1)=X(:,i+1);
        yes=yes+1;
    end
    
end
%plot
figure(1)
plot(t,X(1,:),'r',t,X(2,:),'b',t,X(3,:),'--g',t,y,'k')
xlabel('Time')
ylabel('disp or velo')
grid on

toc