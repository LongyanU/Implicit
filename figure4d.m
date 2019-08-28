clear
clc
close all
% Elapsed time is 33.130840 seconds..
% Elapsed time is 29.353771 seconds.
% new implicit method
nt=403*2;    % number of time steps
isnap=40;    % snapshot sampling

nx=499;
nz=499;

v=ones(nz,nx)*2800;
% v(1:nz/2,:)=1500;

dx=10;  %calculate space increment
h=dx;
x=(0:(nx-1))*dx;
z=(0:(nz-1))*dx;  % initialize space coordinates
dt=0.001; % calculate time step from stability criterion
tau=dt;


f0=65;
t=(1:nt)*dt;
t0=4/f0;                       % initialize time axis
src=10^8*exp(-f0^2*(t-t0).*(t-t0));              % source time function
src=(diff(src))/dx^2;				% time derivative to obtain gaussian


zs=46;
xs=floor(nx/2);
p=zeros([nz nx]); Vx=p; Vz=p;
dz=dx;


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
    
    Px=ifft (ifftshift( M1.*fftshift(fft(p,nx,2),2),2 ),  nx,2);
    Pz=ifft (ifftshift(M2 .*fftshift(fft(p,nz,1),1),1) ,   nz,1);
    Vx=Vx-dt*Px;
    Vz=Vz-dt*Pz;
    [Vx,Vz]=spongeABC(Vx,Vz,nx,nz,45,45,0.009);
    
    
    
    if rem(it,isnap)== 0,
        imagesc(real(p))
        axis equal
        colormap gray
        xlabel('x'),ylabel('z')
        title(sprintf(' Time step: %i - Max ampl: %g ',it,max(max(p))))
        drawnow
    end
    if it==400*2
        pp1=p;
    end
end
toc

save('Figure4dSimulation.mat')