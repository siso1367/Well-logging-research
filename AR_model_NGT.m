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
%~~~~~~~~~~~~~~~~~~~~~~~%
%Main Program for all grade model.
%~~~~~~~~~~~~~~~~~~~~~~~%
zy=[Tho,Ura,Pot,Gamma];
gj=[4,4,5,5];
% g=input('Input the grade model g=');%Input the grade model.
for zz=1:4
    yy=zy(:,zz);
    Y1=yy;
    powerspec1(Y1,10)
    counter=size(Y1);%Characterize size of well logging
    g=gj(:,zz);
    for l=1:g
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
    if(zz==1)
        Q_Th=(theta);%Initial value for the comparison of theta
    else if (zz==3)
            Q_Po=(theta);
        else if (zz==2)
                Q_Ur=(theta);
            end
        end
    end
    %~~~~~~~~~~~~~~~~~~~~~~~%
    %Methods based on least squares.
    %~~~~~~~~~~~~~~~~~~~~~~~%
    n=1:l;
    figure;
    plot(n,AIC,'+');%Akaika criteria.
    title('Akaika criteria');
end
%%
%~~~~~~~~~~~~~~~~~~~~~~~%
%Define a Parametric Model
%~~~~~~~~~~~~~~~~~~~~~~~%
Sampling_time=0.1;
as=1;
%%
%Thorium
bs=Q_Th';
AR_Tho=tf(as,bs,Sampling_time)
[A_Th,B_Th,C_Th,D_Th]=tf2ss(as,bs);
[V_Th,D_Tho] = eig(A_Th);
%%
%Potasium
bs=Q_Po';
AR_Pot=tf(as,bs,Sampling_time)
[A_Po,B_Po,C_Po,D_Po]=tf2ss(as,bs);
[V_Th,D_Pot] = eig(A_Po);
%%
%Uranum
bs=Q_Ur';
AR_Ura=tf(as,bs,Sampling_time);
[A_Ur,B_Ur,C_Ur,D_Ur]=tf2ss(as,bs);
[V_Ur,D_Ura] = eig(A_Ur);
%%
%========================%
%Matlab AR Model
%========================%
Y1=Tho;
m_Th = ar(Y1,4 ,'ls','Ts', Sampling_time)
r_Th = roots(m_Th.a)
Y1=Pot;
m_Po = ar(Y1,5 ,'ls','Ts', Sampling_time)
r_Po = roots(m_Po.a)
Y1=Ura;
m_Ur = ar(Y1,4 ,'ls','Ts', Sampling_time)
r_Ur = roots(m_Ur.a)
Y1=Gamma;
m_Gam = ar(Y1,5 ,'ls','Ts', Sampling_time)
r_Gam = roots(m_Gam.a)
% d = fdesign.notch('N,F0,BW,Ap,Ast',6,.2,.4,.5,60);
% Hd = design(d);
% dsd=Hd.sosMatrix