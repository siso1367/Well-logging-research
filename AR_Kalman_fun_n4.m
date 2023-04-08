% ~~~~~~~~~~~~~~~~~~~~~~~%
% sina soltani 901229.
% ~~~~~~~~~~~~~~~~~~~~~~~%
clear all;
clc;
clf;
close all;
% ~~~~~~~~~~~~~~~~~~~~~~~%
% load system for identification.
% ~~~~~~~~~~~~~~~~~~~~~~~%
Information;
%%
% ~~~~~~~~~~~~~~~~~~~~~~~%
% Matlab AR Model
% ~~~~~~~~~~~~~~~~~~~~~~~%
Sampling_time=0.1;
Y1=gaGamma;
m_Gam = ar(Y1,5 ,'ls','Ts', Sampling_time)
r_Gam = roots(m_Gam.a);

%%
powerspec1(Y1,10)
counter=size(Y1);%Characterize size of well logging
for l=1:10
    for n=1:l %Output
        for k=1:counter
            if(k-n>0)
                z(k,n)=-Y1(k-n,1);%U.
            else
                z(k,n)=0;%U=0;
            end
        end
    end
    theta=[inv(transpose(z)*z)]*[transpose(z)*Y1];%Determine Theta.
    Yo=z*theta;%Output
    e=Y1-Yo;%Error for find grade of model.
    c(n)=det(transpose(z)*z);
    s(n)=transpose(e)*e;
    AIC(n)=1000*log10(s(n))+2*2*n;
end
%~~~~~~~~~~~~~~~~~~~~~~~%
%Methods based on least squares.
%~~~~~~~~~~~~~~~~~~~~~~~%
n=1:l;
figure;
plot(n,AIC,'+');%Akaika criteria.
title('Akaika criteria');

for n=1:5 %Output
        for k=1:counter
            if(k-n>0)
                H(k,n)=Y1(k-n,1);%U.
            else
                H(k,n)=0;%U=0;
            end
        end
end
H=[ones(counter),H];
M=[0,0,0,0,0,0]';
A=eye(6);
Q=2^-20*eye(6);
P = eye(6);
R = [.1];
MM = zeros(size(M,1),size(Y1,2)); PP = zeros(size(M,1),size(M,1),size(Y1,2));
Error_Y=0;
for i=1:counter
    [M,P] = ekf_predict1(M,P,A,Q);
    %%
    Hx(:,:,i)=H(i,:)*eye(6);
    %%
    [M,P,d,s,LH] = ekf_update1(M,P,Y1(i,1),Hx(:,:,i),R);
    MM(:,i) = M;
    PP(:,:,i) = P;
    YY(i)=Hx(:,:,i)*M;
    Err_L2=(Y1(i,1)-YY(i))^2;
    Error_Y=Err_L2+Error_Y;
end
Error_L2=(1/counter(1,1))*(Error_Y)^(1/2);
figure
title('Kalman filtering','FontSize',18);
time=1:counter;
subplot(2,1,1);plot(time,Y1,time,YY,'--r');
subplot(2,1,2);plot(time,MM(3 , :),'r',time,MM(2 , :),'--b',time,MM(1 , :),'-.g',time,MM(4 , :),'--c',time,MM(5 , :),'y',time,MM(6 , :),'--r');
xlabel('Depth','fontsize',12,'color','red');
ylabel('Parameter  %');