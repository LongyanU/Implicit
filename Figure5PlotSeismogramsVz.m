clear;
clc
close all


for iii=195:195
    load('Figure4aSimulation.mat')
    figure;plot(Vz(45:end-45,iii),'b')
    tt1=Vz;
    
    load('Figure4bSimulation.mat')
    hold on;plot(Vz(45:end-45,iii),'m')
    tt2=Vz;
    

    
    load('Figure4dSimulation.mat')
    hold on;plot(real(Vz(45:end-45,iii))   ,'k')
    
    hold on;plot(tt2(45:end-45,iii)-tt1(45:end-45,iii)-0.1,'m')
    legend('Tra implicit', 'Hybrid explicit-implicit','Pseudo spectral method','The difference between the first 2 schemes');
    
    
end
grid on
xlabel('z/dz')
ylabel('Vz(m/s)')
title('')
axis([0 400 -0.52 0.8])
figure;plot(tt2(45:end-45,iii)-tt1(45:end-45,iii),'m')
grid on
xlabel('z/dz')
ylabel('Vz(m/s)')
title('')
axis([0 400 -0.016 0.015])


