clear
clc
close all
% Elapsed time is 33.130840 seconds..
% Elapsed time is 29.353771 seconds.
% new implicit method
nt=403;    % number of time steps
isnap=40;    % snapshot sampling

nx=499;
nz=499;

v=ones(nz,nx)*2800;
% v(1:nz/2,:)=1500;

dx=10;  %calculate space increment
h=dx;
x=(0:(nx-1))*dx;
z=(0:(nz-1))*dx;  % initialize space coordinates
dt=0.002; % calculate time step from stability criterion
tau=dt;


f0=65;
t=(1:nt)*dt;
t0=4/f0;                       % initialize time axis
src=10^8*exp(-f0^2*(t-t0).*(t-t0));              % source time function
src=diff(diff(src))/dx^2;				% time derivative to obtain gaussian


zs=46;
xs=floor(nx/2);

seis_record=zeros(nt,nx);
p=zeros([nz nx]); pnew=p; pold=p;



r=v*dt/h;
coeff=zeros(nz,nx,5);
for ii=1:nz
    for jj=1:nx
        
        coeff(ii,jj,:)=[ 0.507182, 0.156887, -0.00434045, 0.0218515, 0.191447];
        
    end
end


coeff2=cat(3, -coeff(:,:,1)*2  ,coeff(:,:,1)-coeff(:,:,2) ,coeff(:,:,2)-coeff(:,:,3), coeff(:,:,3),coeff(:,:,4) );



b=1-2*coeff(:,:,end);
a = coeff(:,:,end);
cc = a;

for j=1:nx,
    A1{j} = (gallery('tridiag',a(1:nz-1,j),b(:,j),cc(1:nz-1,j)));
end


% a = coeff(end)* ones(nx-1,1);
% cc = a;
A2=cell(nz,1);
for k=1:nz,
    A2 {k}= (gallery('tridiag',a(k,1:nx-1),b(k,:),cc(k,1:nx-1) ));      %%解三对角阵，直接matlab解了
end


Vx=zeros(nz,nx);
Vz=zeros(nz,nx);
d2pzz=p;d2pxx=p;

taper=ones(nz,nx);
for i=1:43
    for j=1:nx
        taper(i,j)=0.5-0.5*cos(pi*(i-1)/(43-1));
        taper(nz-i+1,j)=taper(i,j);
    end
end
for i=1:nz
    for j=1:43
        taper(i,j)=taper(i,j)*(0.5-0.5*cos(pi*(j-1)/(43-1)));
        taper(i,nx-j+1)=taper(i,j);
    end
end

tic
for it=1:nt-2,
    
    
    
    
    d2pz11=(circshift(p,[ 1])+circshift(p,[ -1]));
    d2pz12=(circshift(p,[ 2])+circshift(p,[ -2]));
    d2pz13=(circshift(p,[ 3])+circshift(p,[ -3]));
    
    
    d2px11=(circshift(p,[0  -1])+circshift(p,[0  1]));
    d2px12=(circshift(p,[0  -2])+circshift(p,[0  2]));
    d2px13=(circshift(p,[0  -3])+circshift(p,[0  3]));
    
    
    commonPoints=circshift(p,[-1 -1])+circshift(p,[1  1]) +circshift(p,[1 -1])+circshift(p,[-1 1]);
    d2pzAddedpoints= -2*d2px11 +commonPoints;
    d2pxAddedpoints=-2*d2pz11+commonPoints;
    
    d2px=coeff2(:,:,1).*p+coeff2(:,:,2).*d2px11+coeff2(:,:,3).*d2px12+coeff2(:,:,4).*d2px13+coeff2(:,:,5).*d2pxAddedpoints;
    d2pz=coeff2(:,:,1).*p+coeff2(:,:,2).*d2pz11+coeff2(:,:,3).*d2pz12+coeff2(:,:,4).*d2pz13 +coeff2(:,:,5).*d2pzAddedpoints;
    
    for j=1:nx,
        d2pzz(:,j)=A1{j}\d2pz(:,j);   %%解三对角阵，直接matlab解了
    end
    for k=1:nz,
        d2pxx(k,:)=A2{k}\d2px(k,:)';       %%解三对角阵，直接matlab解了
    end
    
    pnew=2*p-pold+v.*v.*(d2pzz +d2pxx)*dt^2/dx^2;
    
    for kkk=1:20
        pdan=zeros(size(p));
        %底边
        k=nz-30+kkk;
        for j=32+1-kkk:nx-31+kkk-1
            pboundarynew(k,j)= (2*h*tau^2*v(k,j)*((pold(k,j) - pold(k-1,j) +  pnew(k-1,j))/(2*h*tau) - (pold(k,j) - 2*p(k-1,j)- 2*p(k,j) + pold(k-1,j) ...
                + pnew(k-1,j))/(2*tau^2*v(k,j)) + (v(k,j)*(pold(k,j+1) - 2*pold(k,j) + pold(k,j-1) - 2*pnew(k-1,j) + pnew(k-1,j+1)+ pnew(k-1,j-1)))/(4*h^2)))/(h + tau*v(k,j));
        end
        
        % 左边界
        j=30-kkk+1;
        for k=32+1-kkk:nz-31+kkk-1
            pboundarynew(k,j)=(2*h*tau^2*v(k,j)*((pold(k,j) - pold(k,j+1) + pnew(k,j+1))/(2*h*tau) - (pold(k,j) - 2*p(k,j+1) - 2*p(k,j) + pold(k,j+1) + pnew(k,j+1))/(2*tau^2*v(k,j)) +...
                (v(k,j)*(pold(k+1,j) - 2*pold(k,j) + pold(k-1,j) - 2*pnew(k,j+1) + pnew(k+1,j+1) + pnew(k-1,j+1)))/(4*h^2)))/(h + tau*v(k,j));
        end
        
        %右边界
        j=nx-30+kkk;
        for k=32+1-kkk:nz-31+kkk-1
            pboundarynew(k,j)=(2*h*tau^2*v(k,j)*((pold(k,j) - pold(k,j-1) + pnew(k,j-1))/(2*h*tau) - (pold(k,j) - 2*p(k,j-1) - 2*p(k,j) + pold(k,j-1) + pnew(k,j-1))/(2*tau^2*v(k,j)) ...
                + (v(k,j)*(pold(k+1,j) - 2*pold(k,j) + pold(k-1,j) - 2*pnew(k,j-1) + pnew(k+1,j-1) + pnew(k-1,j-1)))/(4*h^2)))/(h + tau*v(k,j));
        end
        %
        
        k=30-kkk+1;
        for j=32+1-kkk:nx-31+kkk-1
            pboundarynew(k,j)=(2*h*tau^2*v(k,j)*((pold(k,j) - pold(k+1,j) + pnew(k+1,j))/(2*h*tau) - (pold(k,j) - 2*p(k+1,j) - 2*p(k,j) + pold(k+1,j) + pnew(k+1,j))/(2*tau^2*v(k,j)) + ...
                (v(k,j)*(pold(k,j+1) - 2*pold(k,j) + pold(k,j-1) - 2*pnew(k+1,j) + pnew(k+1,j+1) + pnew(k+1,j-1)))/(4*h^2)))/(h + tau*v(k,j));
        end
        
        
        %右上角的三个边界点
        k=30-kkk+1;
        j=nx-30+kkk;
        %j-1
        pboundarynew(k,j-1)=pnew(k+1,j-1)/h+pnew(k,j-1-1)/h+p(k,j-1)*sqrt(2)/v(k,j)/tau;
        pboundarynew(k,j-1)=pboundarynew(k,j-1)/(2/h+sqrt(2)/v(k,j)/tau);
        %k-1
        pboundarynew(k+1,j)=pnew(k+1+1,j)/h+pnew(k+1,j-1)/h+p(k+1,j)*sqrt(2)/v(k,j)/tau;
        pboundarynew(k+1,j)=pboundarynew(k+1,j)/(2/h+sqrt(2)/v(k,j)/tau);
        
        pboundarynew(k,j)=pnew(k+1,j)/h+pnew(k,j-1)/h+p(k,j)*sqrt(2)/v(k,j)/tau;
        pboundarynew(k,j)=pboundarynew(k,j)/(2/h+sqrt(2)/v(k,j)/tau);
        
        
        %左上角的三个边界点
        k=30-kkk+1;
        j=30-kkk+1;
        %j-1
        pboundarynew(k,j+1)=pnew(k+1,j+1)/h+pnew(k,j+1+1)/h+p(k,j+1)*sqrt(2)/v(k,j)/tau;
        pboundarynew(k,j+1)=pboundarynew(k,j+1)/(2/h+sqrt(2)/v(k,j)/tau);
        
        %k-1
        pboundarynew(k+1,j)=pnew(k+1+1,j)/h+pnew(k+1,j+1)/h+p(k+1,j)*sqrt(2)/v(k,j)/tau;
        pboundarynew(k+1,j)=pboundarynew(k+1,j)/(2/h+sqrt(2)/v(k,j)/tau);
        
        pboundarynew(k,j)=pnew(k+1,j)/h+pnew(k,j+1)/h+p(k,j)*sqrt(2)/v(k,j)/tau;
        pboundarynew(k,j)=pboundarynew(k,j)/(2/h+sqrt(2)/v(k,j)/tau);
        
        
        %左下角 这个还是很有作用的
        k=nz-30+kkk;
        j=30-kkk+1;
        
        %j+1
        pboundarynew(k,j+1)=pnew(k-1,j+1)/h+pnew(k,j+1+1)/h+p(k,j+1)*sqrt(2)/v(k,j)/tau;
        pboundarynew(k,j+1)=pboundarynew(k,j+1)/(2/h+sqrt(2)/v(k,j)/tau);
        
        %k-1
        pboundarynew(k-1,j)=pnew(k-1-1,j)/h+pnew(k-1,j+1)/h+p(k-1,j)*sqrt(2)/v(k,j)/tau;
        pboundarynew(k-1,j)=pboundarynew(k-1,j)/(2/h+sqrt(2)/v(k,j)/tau);
        
        pboundarynew(k,j)=pnew(k-1,j)/h+pnew(k,j+1)/h+p(k,j)*sqrt(2)/v(k,j)/tau;
        pboundarynew(k,j)=pboundarynew(k,j)/(2/h+sqrt(2)/v(k,j)/tau);
        
        k=nz-30+kkk;
        j=nx-30+kkk;
        %j-1
        pboundarynew(k,j-1)=pnew(k-1,j-1)/h+pnew(k,j-1-1)/h+p(k,j-1)*sqrt(2)/v(k,j)/tau;
        pboundarynew(k,j-1)=pboundarynew(k,j-1)/(2/h+sqrt(2)/v(k,j)/tau);
        
        %k-1
        pboundarynew(k-1,j)=pnew(k-1-1,j)/h+pnew(k-1,j-1)/h+p(k-1,j)*sqrt(2)/v(k,j)/tau;
        pboundarynew(k-1,j)=pboundarynew(k-1,j)/(2/h+sqrt(2)/v(k,j)/tau);
        
        pboundarynew(k,j)=pnew(k-1,j)/h+pnew(k,j-1)/h+p(k,j)*sqrt(2)/v(k,j)/tau;
        pnew(k,j)=0.1*kkk*pdan(k,j)+(1-0.1*kkk)*pnew(k,j);
        %右下角
    end
    
    for kkk=1:20
        pnew(nz-30+kkk,30+1-kkk:nx-29+kkk-1)=0.05*kkk*pboundarynew(nz-30+kkk,30+1-kkk:nx-29+kkk-1)+(1-0.05*kkk)*pnew(nz-30+kkk,30+1-kkk:nx-29+kkk-1);
        pnew(30+1-kkk:nz-29+kkk-1,30-kkk+1)=0.05*kkk*pboundarynew(30+1-kkk:nz-29+kkk-1,30-kkk+1)+(1-0.05*kkk)*pnew(30+1-kkk:nz-29+kkk-1,30-kkk+1);
        pnew(30+1-kkk:nz-29+kkk-1,nx-30+kkk)=0.05*kkk*pboundarynew(30+1-kkk:nz-29+kkk-1,nx-30+kkk)+(1-0.05*kkk)*pnew(30+1-kkk:nz-29+kkk-1,nx-30+kkk);
        pnew(30-kkk+1,30+1-kkk:nx-29+kkk-1)=0.05*kkk*pboundarynew(30-kkk+1,30+1-kkk:nx-29+kkk-1)+(1-0.05*kkk)*pnew(30-kkk+1,30+1-kkk:nx-29+kkk-1);
    end
    
    
    pnew(zs,xs)=pnew(zs,xs)+src(it);
    pold=p;											% time lev(k,j)els
    p=pnew;
    
    [p,pold]=spongeABC(p,pold,nx,nz,45,45,0.0009);
    
    if rem(it,isnap)== 0,
        it
        %         imagesc(x,z,p,[-0.0007 0.0007]), axis equal
        imagesc(x,z,p), axis equal
        colormap gray
        xlabel('x'),ylabel('z')
        title(sprintf(' Time step: %i - Max ampl: %g ',it,max(max(p))))
        drawnow
    end
    
    if it==400
        pp1=p;
    end
    
end


toc
save('Figure4cSimulation.mat')
figure;imagesc(pp1(45:end-45,45:end-45))
xlabel('x/dx')
ylabel('z/dz')
h = colorbar;