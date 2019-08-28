clear;
clc
close all;

load('Figure7aSaltTraScheme.mat')
figure;plot([1:2500]*2,-real(seis_recordVz(1:1:2500,450)),'b')
tt1=-seis_recordVz(1:2500,450);

load('Figure7bSaltNewScheme.mat')
hold on ;plot([1:2500]*2,-real(seis_recordVz(1:1:2500,450)),'k')
tt2=-seis_recordVz(1:2500,450);


load('Figure7dKSpace.mat')
hold on;plot([1:2500]*2,-real(seis_recordVz(1:2:5000,450)),'r')


hold on ;plot([1:2500]*2,tt2-tt1-0.05,'k')

axis([0 5000 -0.4 0.5])
grid on
legend('Tra implicit scheme','Explicit-implicit scheme', 'Pseudo-spectrum method','the difference bewteen the first 2 schemes')
ylabel('Vz(m/s)')
xlabel('travel time(ms)')

