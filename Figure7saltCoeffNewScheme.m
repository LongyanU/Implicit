clear;
clc
close all

global v ratio M h tau
h=15;
tau=0.002;
M=3;
% Elapsed time is 989.375542 seconds.

x0=0.001*ones(1,M+2);%系数的初值是0

options = optimset('Algorithm','levenberg-marquardt','TolFun',10^-10,'TolX',10^-10,'MaxFunEvals',2000,'MaxIter',200);


coeffJune=zeros(4790-1486+1,7);
ii=1;
tic
for v=1486:4790
    v

        ratio=0.82;

        x0=0.001*ones(1,M+2);
        aa=-5*ones(1,M+1);
        bb=5*ones(1,M+1);
        [x,resnorm] = lsqnonlin(@myfun7,x0,aa,bb,options) ;
        April19(ii,:)=real(x);
        %     else
        %         ratio=0.5;
        %        [x,resnorm]     = lsqnonlin(@myfun7,x0,aa,bb,options) ;   % Invoke optimizer

    ii=ii+1;
end
toc
save('coeffApril19.mat','April19')