clear;
clc
close all


for iii=195:195
    load('Figure4aSimulation.mat')
    figure;plot(Vx(45:end-45,iii),'b')
    tt1=Vx;
    
    load('Figure4bSimulation.mat')
    hold on;plot(Vx(45:end-45,iii),'m')
    tt2=Vx;
    

    
    load('Figure4dSimulation.mat')
    hold on;plot(real(Vx(45:end-45,iii))   ,'k')
    
    hold on;plot(tt2(45:end-45,iii)-tt1(45:end-45,iii)-2*10^-2,'m')
    legend('Tra implicit', 'Hybrid explicit-implicit','Pseudo spectral method','The difference between the first 2 schemes');
    
    
end
grid on
xlabel('z/dz')
ylabel('Vx(m/s)')
title('')
axis([0 400 -0.15 0.2])
figure;plot(tt2(45:end-45,iii)-tt1(45:end-45,iii),'m')
grid on
xlabel('z/dz')
ylabel('Vx(m/s)')
title('')
% axis([0 400 -1.2*10^-4 1*10^-4])


