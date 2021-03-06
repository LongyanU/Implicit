clear;
clc
close all;


load('Figure7aSaltTraScheme.mat')
plotimage(-seis_recordp(:,45:end-45))
xlabel('x/dx')
ylabel('Travel time(ms)')
title('')


load('Figure7bSaltNewScheme.mat')
plotimage(-seis_recordp(:,45:end-45))
xlabel('x/dx')
ylabel('Travel time(ms)')
title('')


load('Figure7cSecondOrderScheme.mat')
plotimage(-seis_recordp(:,45:end-45))
xlabel('x/dx')
ylabel('Travel time(ms)')
title('')


load('Figure7dKSpace.mat')
plotimage(-real(seis_recordp(:,45:end-45)))
xlabel('x/dx')
ylabel('Travel time(ms)')
title('')


