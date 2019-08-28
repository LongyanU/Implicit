clear;
clc
close all;

load('Figure7aSaltTraScheme.mat')
figure;plot([1:2500]*2,-real(seis_recordVx(1:1:2500,450)),'b')
tt1=-seis_recordVx(1:2500,450);

load('Figure7bSaltNewScheme.mat')
hold on ;plot([1:2500]*2,-real(seis_recordVx(1:1:2500,450)),'k')
tt2=-seis_recordVx(1:2500,450);


load('Figure7dKSpace.mat')
hold on;plot([1:2500]*2,-real(seis_recordVx(1:2:5000,450)),'r')


hold on ;plot([1:2500]*2,tt2-tt1-0.1,'k')

axis([0 5000 -5.2 3.5])
grid on
legend('Tra implicit scheme','Explicit-implicit scheme', 'Pseudo-spectrum method','the difference bewteen the first 2 schemes')
ylabel('Vx(m/s)')
xlabel('travel time(ms)')
hold on ;plot([1:2500]*2,tt2-tt1-0.015,'k')
hold on ;plot([1:2500]*2,tt2-tt1-0.003,'k')

