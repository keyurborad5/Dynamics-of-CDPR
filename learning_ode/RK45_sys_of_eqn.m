clc
clear all
%defining input force
time=0:0.1:30;
fmea=5*rectangularPulse(2,5,time)+10*rectangularPulse(10,14,time)+3*rectangularPulse(20,22,time);
%fmea=10;
fd=0;

%% define functoin handle
%X=[disp;velo] ==> disp=X(1,:),velo=X(2,:)
%i made function handle with so many variables because i need to change its
%value while in a loop of solving RK45
f=@(t,X,k,m,b,x_ref,fmea)[[0 1;-k/m -b/m ]*[X(1);X(2)]+[0;(fmea-fd+k*x_ref)/m]];

%% initial conditions
t(1)=0;
X(:,1)=[0,0];
%% step size
h=0.1;
tfinal=30;
N=tfinal/h;
x_ref=0;
%% RK45
for i=1:N
    m=3;b=4;
    %This conditional statements deals with being stiff and compliant
    if X(1,i)< 3+x_ref
        k=2.5;%behave stiff
        
    else
        k=0;%behave compliant
        if X(1,i)-X(1,i-1)<1e-2 %checks if the position is constant for some time
            count=count+1;
            if count >=20 % that some time is decided by this no. of counts
                x_ref=X(1,i); %changes the reference around which it needed to be stiff
            end
        else
            count=0;
        end
    end
    %update time
    t(i+1)=t(i)+h;
    %upsade of X
    k1=f(t(i)    ,X(:,i)       ,k,m,b,x_ref,fmea(i));
    k2=f(t(i)+h/2,X(:,i)+k1*h/2,k,m,b,x_ref,fmea(i));
    k3=f(t(i)+h/2,X(:,i)+k2*h/2,k,m,b,x_ref,fmea(i));
    k4=f(t(i)+h  ,X(:,i)+k3*h  ,k,m,b,x_ref,fmea(i));
    X(:,i+1)=X(:,i)+(h/6)*(k1+2*k2+2*k3+k4);
end
%plot
figure(1)
plot(t,X(1,:),'r',t,X(2,:),'--b',t,fmea,'k')
xlabel('Time')
ylabel('disp or velo')
grid on
%legend('disp','velocity')