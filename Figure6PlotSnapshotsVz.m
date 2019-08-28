 clear;
 clc
 close all


load('Figure4aSimulation.mat')
plotimage(Vx(45:end-45,45:end-45))
xlabel('x/dx')
ylabel('z/dz')
colormap gray
tttt1=Vz;
h = colorbar;

load('Figure4bSimulation.mat')
plotimage((real(Vx(45:end-45,45:end-45))))
xlabel('x/dx')
ylabel('z/dz')
colormap gray
h = colorbar;

load('Figure4dSimulation.mat')
plotimage((real(Vx(45:end-45,45:end-45))))
xlabel('x/dx')
ylabel('z/dz')
h = colorbar;
colormap gray
h = colorbar;