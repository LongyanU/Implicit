 clear;
 clc
 close all


load('Figure4aSimulation.mat')
plotimage(p(45:end-45,45:end-45))
xlabel('x/dx')
ylabel('z/dz')
h = colorbar;
colormap gray
tttt1=Vz;

load('Figure4bSimulation.mat')
plotimage((real(p(45:end-45,45:end-45))))
xlabel('x/dx')
ylabel('z/dz')
h = colorbar;
colormap gray

load('Figure4cSimulation.mat')
plotimage((real(p(45:end-45,45:end-45))))
xlabel('x/dx')
ylabel('z/dz')
h = colorbar;
colormap gray

load('Figure4dSimulation.mat')
plotimage((real(p(45:end-45,45:end-45))))
xlabel('x/dx')
ylabel('z/dz')
h = colorbar;
colormap gray