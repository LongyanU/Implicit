% first run Figure7aSaltTraScheme.m,Figure7bSaltNewScheme.m,figure7cSecondOrderScheme.m,figure7dKSpace.m to do the
% wave equation simulation, 
% then run Figure8CompareSeismicRecordsP.m, figure8CompareSeismogramsP.m, figure9CompareSeismicRecordsVx.m,
% figure9CompareSeismogramsVx.m,figure10CompareSeismicRecordsVz.m,figure10CompareSeismogramsVz.m to get the seismic records.
clear
clc %%%%%%%
close all
nt=2502*2;    % number of time steps
eps=.6;     % stability
isnap=50;    % snapshot sampling
load('vv')

c1=flipud(c);

v=c1;
nx=799;
nx=nx+45*2;
nz=475;
nz=nz+45*2;

vv=zeros(nz,nx);
for ii=1:nz-90
    for jj=1:nx-90
        vv(ii+45,jj+45)=v(ii,jj);
    end
end

for ii=1:nz-90  %%left
    for jj=1:45
        vv(ii+45,jj)=v(ii,1);
    end
end

for ii=1:nz-90  %%right
    for jj=nx-45:nx
        vv(ii+45,jj)=v(ii,800);
    end
end


for ii=1:45  %%top
    for jj=1:nx
        vv(ii,jj)=vv(46,jj);
    end
end

for ii=nz-44:nz  %%bottom
    for jj=1:nx
        vv(ii,jj)=vv(nz-45,jj);
    end
end


clear v
v=vv;
% v=round(v);


dx=15;  %calculate space increment
h=dx;
dz=dx;
x=(0:(nx-1))*dx;
z=(0:(nz-1))*dx;  % initialize space coordinates
dt=0.001; % calculate time step from stability criterion
tau=dt;


f0=45;
t=(1:nt)*dt;
t0=4/f0;                       % initialize time axis
src=10^8*exp(-f0^2*(t-t0).*(t-t0));              % source time function
src=(diff(src))/dx^2;				% time derivative to obtain gaussian


zs=46;
xs=600-150+25;

seis_recordVx=zeros(nt,nx);
seis_recordVz=zeros(nt,nx);
seis_recordp=zeros(nt,nx);
p=zeros([nz nx]); pnew=p; pold=p;

VxzAddedpoint=p;

r=v*dt/h;

Vx=zeros(nz,nx);
Vz=zeros(nz,nx);
d2pzz=p;d2pxx=p;

kx=linspace(-pi/dx,pi/dx, nx);
kz=linspace(-pi/dz,pi/dz, nz);

kexp=1i*kx.*exp(1i*kx*dx/2);
kexpp=(1i*kz.*exp(1i*kz*dx/2))';
kexpm= 1i*kx.*exp(-1i*kx*dx/2);
kexppm=(1i*(kz.*exp(-1i*kz*dx/2)))';
M1=repmat(kexp,nz,1);
M2=repmat(kexpp,1,nx);
M3=repmat(kexpm,nz,1);
M4=repmat(kexppm,1,nx);
tic
for it=1:nt-2
    
    vxx=ifft( ifftshift(M3.*fftshift(fft(Vx,nx,2),2 ),2) ,nx,2);
    vzz=ifft(ifftshift(M4.*fftshift(fft(Vz,nz,1),1),1), nz,1);
    p=p-dt*v.^2.*(vzz+vxx);
    
    [p,p]=spongeABC(p,p,nx,nz,45,45,0.009);
    p(zs,xs)=p(zs,xs)+src(it);
    seis_recordp(it,:)=p(zs,:);
    
    Px=ifft (ifftshift( M1.*fftshift(fft(p,nx,2),2),2 ),  nx,2);
    Pz=ifft (ifftshift(M2 .*fftshift(fft(p,nz,1),1),1) ,   nz,1);
    Vx=Vx-dt*Px;
    Vz=Vz-dt*Pz;
    [Vx,Vz]=spongeABC(Vx,Vz,nx,nz,45,45,0.009);
    
    seis_recordVx(it,:)=Vx(zs,:);
    seis_recordVz(it,:)=Vz(zs,:);
    
    if rem(it,isnap)== 0,
        imagesc(real(p))
        axis equal
        colormap gray
        xlabel('x'),ylabel('z')
        title(sprintf(' Time step: %i - Max ampl: %g ',it,max(max(p))))
        drawnow
    end
    
end
toc

save Figure7dKSpace.mat