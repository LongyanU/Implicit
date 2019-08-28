clear;
clc
close all;


load('Figure7aSaltTraScheme.mat')
figure;plot(-seis_recordp(1:end,450),'b')
tt1=-seis_recordp(1:end,450);

load('Figure7bSaltNewScheme.mat')
hold on ;plot(-seis_recordp(1:end,450),'c')
tt2=-seis_recordp(1:end,450);
hold on ;plot(tt2-tt1-400,'k')


load('Figure7cSecondOrderScheme.mat')
hold on;plot(-seis_recordp(1:end,450),'m')

load('Figure7dKSpace.mat')
hold on;plot(-real(seis_recordp(1:2:end,450)),'r')


axis([0 2450 -4800 7600])
grid on
legend('Tra implicit scheme','Explicit-implicit scheme','the difference bewteen the first 2 schemes', 'Second-order implicit scheme', 'Pseudo-spectrum method')
ylabel('p (Pa)')
xlabel('travel time(m/s)')
hold on ;plot(tt2-tt1-30,'k')
