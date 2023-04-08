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
fs=60;
counter=size(X1);
%%
%~~~~~~~~~~~~~~~~~~~~~~~%
%Normalize Thorium Data.
%~~~~~~~~~~~~~~~~~~~~~~~%
tho=Tho;
nor=1/counter(1,1):1/counter(1,1):1;
wavenum=max(tho);
for i=1:counter(1,1)
    if tho(i)==wavenum
        num_tho=i/counter(1,1)
    end
end
% M=mean(tho);
% tho=tho-M;
% hold on
% cs = csapi(nor,tho);
% fnplt(cs);
% hold on;
% sp = spap2(4,5,nor,tho);
% fnplt(sp,'g');
figure
spline_ls=LS_spline(nor,tho);
Tho_spline_ls=spline_ls;
m=max(Tho_spline_ls);
Tho_spline_ls=Tho_spline_ls/m;
power000=powerspec2(Tho_spline_ls,fs);
figure;
plot(nor,Tho_spline_ls,'r');
hold on
x=Tho_spline_ls;
a=.4;
y(1)=x(1);y(2)=x(2);
for ii=1:100
    for i=3:counter(1,1)
        y(i)=x(i)-2*cos(pi*num_tho)*x(i-1)+x(i-2)+2*a*cos(pi*num_tho)*y(i-1)-a^2*y(i-2);
    end
    wavenum=max(y);
    for i=1:counter(1,1)
        if tho(i)==wavenum
            num_tho=i/counter(1,1)
        end
    end
end
plot(nor,y)
%%
%~~~~~~~~~~~~~~~~~~~~~~~%
%Normalize Uranium Data.
%~~~~~~~~~~~~~~~~~~~~~~~%
ura=Ura;
nor=1/counter(1,1):1/counter(1,1):1;
wavenum=max(ura);
for i=1:counter(1,1)
    if ura(i)==wavenum
        num_ura=i/counter(1,1)
    end
end
% M=mean(ura);
% ura=ura-M;
figure
spline_ls=LS_spline(nor,ura);
Ur_spline_ls=spline_ls;
m=max(Ur_spline_ls);
Ur_spline_ls=Ur_spline_ls/m;
figure;
plot(nor,Ur_spline_ls,'r');
hold on
x=Ur_spline_ls;
a=.4;
y(1)=x(1);y(2)=x(2);
for ii=1:100
    for i=3:counter(1,1)
        y(i)=x(i)-2*cos(pi*num_tho)*x(i-1)+x(i-2)+2*a*cos(pi*num_tho)*y(i-1)-a^2*y(i-2);
    end
    wavenum=max(y);
    for i=1:counter(1,1)
        if tho(i)==wavenum
            num_tho=i/counter(1,1)
        end
    end
end
plot(nor,y)
powerspec2(Ur_spline_ls,fs);
%%
%~~~~~~~~~~~~~~~~~~~~~~~%
%Normalize Potassium Data.
%~~~~~~~~~~~~~~~~~~~~~~~%
pot=Pot;
nor=1/counter(1,1):1/counter(1,1):1;
wavenum=max(pot);
for i=1:counter(1,1)
    if pot(i)==wavenum
        num_pot=i/counter(1,1)
    end
end
% M=mean(pot);
% pot=pot-M;
figure;
spline_ls=LS_spline(nor,pot);
Pot_spline_ls=spline_ls;
m=max(Pot_spline_ls);
Pot_spline_ls=Pot_spline_ls/m;
figure;
plot(nor,Pot_spline_ls,'r');
hold on
x=Pot_spline_ls;
a=.4;
y(1)=x(1);y(2)=x(2);
for ii=1:100
    for i=3:counter(1,1)
        y(i)=x(i)-2*cos(pi*num_tho)*x(i-1)+x(i-2)+2*a*cos(pi*num_tho)*y(i-1)-a^2*y(i-2);
    end
    wavenum=max(y);
    for i=1:counter(1,1)
        if tho(i)==wavenum
            num_tho=i/counter(1,1)
        end
    end
end
plot(nor,y)
powerspec2(Pot_spline_ls,fs);
%%
%~~~~~~~~~~~~~~~~~~~~~~~%
%Normalize Gamma Data.
%~~~~~~~~~~~~~~~~~~~~~~~%
gamma=Gamma;
nor=1/counter(1,1):1/counter(1,1):1;
wavenum=max(gamma);
for i=1:counter(1,1)
    if gamma(i)==wavenum
        num_gamma=i/counter(1,1)
    end
end
% M=mean(gamma);
% gamma=gamma-M;
% figure;
% plot(nor,gamma);
figure;
spline_ls=LS_spline(nor,gamma);
Gamma_spline_ls=spline_ls;
m=max(Gamma_spline_ls);
Gamma_spline_ls=Gamma_spline_ls/m;
figure;
plot(nor,Gamma_spline_ls,'r');
hold on
x=Gamma_spline_ls;
a=.4;
y(1)=x(1);y(2)=x(2);
for ii=1:100
    for i=3:counter(1,1)
        y(i)=x(i)-2*cos(pi*num_tho)*x(i-1)+x(i-2)+2*a*cos(pi*num_tho)*y(i-1)-a^2*y(i-2);
    end
    wavenum=max(y);
    for i=1:counter(1,1)
        if tho(i)==wavenum
            num_tho=i/counter(1,1)
        end
    end
end
plot(nor,y)
powerspec2(Gamma_spline_ls,fs);
%%
%~~~~~~~~~~~~~~~~~~~~~~~%
%Normalize Caliper Data.
%~~~~~~~~~~~~~~~~~~~~~~~%
caliper=Caliper;
nor=1/counter(1,1):1/counter(1,1):1;
% M=mean(caliper);
% caliper=caliper-M;
wavenum=max(caliper);
for i=1:counter(1,1)
    if caliper(i)==wavenum
        num_caliper=i/counter(1,1)
    end
end
% figure;
% hold on
% cs = csapi(nor,caliper);
% fnplt(cs);
% hold on;
% sp = spap2(3,4,nor,caliper);
% fnplt(sp,'g.');
% hold on;
% plot(nor,caliper);
% % cs=cs-sp;
% % hold on;
% % plot(nor,caliper);
% powerspec(tho,fs);
figure;
spline_ls=LS_spline(nor,caliper);
Caliper_spline_ls=spline_ls;
m=max(Caliper_spline_ls);
Caliper_spline_ls=Caliper_spline_ls/m;
figure;
plot(nor,Caliper_spline_ls,'r');
hold on
x=Caliper_spline_ls;
a=.4;
y(1)=x(1);y(2)=x(2);
for ii=1:100
    for i=3:counter(1,1)
        y(i)=x(i)-2*cos(pi*num_tho)*x(i-1)+x(i-2)+2*a*cos(pi*num_tho)*y(i-1)-a^2*y(i-2);
    end
    wavenum=max(y);
    for i=1:counter(1,1)
        if tho(i)==wavenum
            num_tho=i/counter(1,1)
        end
    end
end
plot(nor,y)
powerspec2(Caliper_spline_ls,fs);
powerspec2(Caliper_spline_ls,fs);
%%
%~~~~~~~~~~~~~~~~~~~~~~~%
%subtract Data.
%~~~~~~~~~~~~~~~~~~~~~~~%
% splinetool(nor,caliper)