clear;
clc
close all


for iii=195:195
    load('Figure4aSimulation.mat')
    figure;plot(p(45:end-45,iii),'b')
    tt1=p;
    
    load('Figure4bSimulation.mat')
    hold on;plot(p(45:end-45,iii),'m')
    tt2=p;
    
    load('Figure4cSimulation.mat')
    hold on;plot(real(p(45:end-45,iii))   ,'r')
    
    load('Figure4dSimulation.mat')
    hold on;plot(real(p(45:end-45,iii))   ,'k')
    
    hold on;plot(tt2(45:end-45,iii)-tt1(45:end-45,iii)-500,'m')
    legend('Tra implicit', 'Hybrid explicit-implicit','Implicit scheme for second-order wave equation','Pseudo spectral method','The difference between the first 2 schemes');
    
    
end
grid on
xlabel('z/dz')
ylabel('p (Pa)')
title('')
axis([0 400 -2300 1500])
figure;plot(tt2(45:end-45,iii)-tt1(45:end-45,iii),'m')
grid on
xlabel('z/dz')
ylabel('p (Pa)')
title('')
axis([0 400 -32 30])


