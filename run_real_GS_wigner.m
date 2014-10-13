
%% parameter settings,
wo=9; % reference beam waist 
A=1; % amplitude constant 
z=100; % lateral distance from the beam waist
lambda=500e-6; % wavelength in millimeters
lc=7; % lateral coherent length in millimeters

nx=512;   % sum of sampling points in x direction
max_x=20; % maximum value on x axis
nu=nx;    % sum of sampling points in u direction
max_u=6.375;%1;%0.1; % maximum value on u axis, but the maximum of u should be lower than 1/dx/2 according to the sampling theory of wigner distribution


%% generate the gaussian schell beam in wigner space
% real units in millimeters are used 
[GS_Corr,Wig_GS]=func_gnr_GS_wig(nx,max_x,nu,max_u,wo,z,lambda,lc);



%% read in the Fresnel lens from Gashaw's data

zF=importdata('D:\ILLUMINATION PROJECT\Matlab-Wigner\Fresnel lens from Diinesh\Gashaw\zF123p.mat'); %zF123p
zFL=zF(1:(length(zF)-2));

yF=importdata('D:\ILLUMINATION PROJECT\Matlab-Wigner\Fresnel lens from Diinesh\Gashaw\yF123p.mat'); %yF123p

yFL=yF(1:(length(zF)-2));

max_FL=max(yFL);

n=nx/2;

yyFL=linspace(0,max_FL,n);

zzFL=interp1(yFL,zFL,yyFL);% interp1, spline

zzFL(isnan(zzFL))=zF(1,1);

zzFL=zzFL-min(zzFL);


Fi=zeros(1,2*n);
Fi(1,(n+1):(n*2))=zzFL;
Fi(1,1:n)=fliplr(zzFL);

FL=-Fi+1+2.375;

figure(6)
plot(FL,'.')
title('Fresnel lens')
% axis([min(xi),max(xi),-5,5]);

%% generate the 2D grid of x1, x2 for futher calculation of cross-corelation function

[xi,ui]=gnr_array(nx,nx,max_x,max_u);

[x,xx] = meshgrid(xi, xi);%(xi,xxi);
% % x changes in columns, xx changes in lines
[xx2,u]=meshgrid(xi,ui);%(xxi,ui);
% %xx2 is actually xx, but changes in columns. u is the frequencies changing in
% %lines.

x1=x+xx/2;
x2=x-xx/2;


x1(x1>max(xi))=0; % it gives the same result as F1(isnan(F1))=0;
x1(x1<min(xi))=0;

% two functions for cross-correlation
F1=interp1(xi,FL,x1);
F1(isnan(F1))=0;

F2=interp1(xi,FL,x2);
F2(isnan(F2))=0;
%% wigner calculation

dx=abs(xi(2)-xi(1));

n=1.5;
lambda=550e-6;
F_phase1=exp(1i*2*pi*(n-1)/lambda.*F1);
F_phase2=exp(1i*2*pi*(n-1)/lambda.*F2);

% Corr=GS_Corr.*F_phase1.*conj(F_phase2);
Corr=F_phase1.*conj(F_phase2);
C=xcorr2(F_phase1,F_phase2);
% Corr=GS_Corr;


Phase=exp(-2i*pi*u.*xx2);
% corresponding-pixel multiplication 
% phase shift factor, horizontal axis is xx (rr), vertical axis is q.

Wig=Phase*Corr*dx;
%%
figure(3)
pcolor(xi,ui,real(Wig));
colormap jet, shading interp
colorbar('location','eastoutside')
set(gca,'FontSize',14,'FontSize',20);
title('Wig(x,u)','FontSize',20)
xlabel('x','FontSize',20)
ylabel('u','FontSize',20)

%
figure(4)
pcolor(xi,xi,real(GS_Corr));
colormap jet, shading interp
colorbar('location','eastoutside')
set(gca,'FontSize',14,'FontSize',20);
title('GauSchell Corr')
xlabel('x or x_1-x_2','FontSize',20)
ylabel('xx or (x_1+x_2)/2)','FontSize',20)

figure(5)
pcolor(imag(F_phase1.*conj(F_phase2)));
colormap jet, shading interp
colorbar('location','eastoutside')
set(gca,'FontSize',14,'FontSize',20);
title('real(FLphase(x1).*conj(FLphase(x2))')
xlabel('x or x_1-x_2','FontSize',20)
ylabel('xx or (x_1+x_2)/2)','FontSize',20)
%
Wig_conv=conv2(Wig_GS,Wig);

figure(7)
pcolor(real(Wig_conv));
colormap jet, shading interp
colorbar('location','eastoutside')
set(gca,'FontSize',14,'FontSize',20);
title('Wig conv')
xlabel('x','FontSize',20)
ylabel('u','FontSize',20)

plot_2d(ui,'u', 

