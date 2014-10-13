function [Corr_rect_GS,Wig_rect]=func_gnr_rectGS_wig(nx,max_x,lambda,M,r,sigma,delta,theta,path)

% M=40;

nu=nx;      % sum of sampling points in u direction
dx = max_x*2 / (nx-2);
max_u=1/dx/2*lambda/2; %0.0318  unit: radian, 1/dx/2: half of a fourier transform, scale: lambda/2: half of a period

% [xi,ui]=gnr_array(nx,nu,max_x,max_u);

xi=gnr_array(nx,max_x);

ui=gnr_array(nu,max_u);

x_prime=xi*2;

[x,xx] = meshgrid(xi, x_prime);%(xi,xxi);
% % x changes in columns, xx changes in lines
[xx2,u]=meshgrid(x_prime,ui);%(xxi,ui);
% %xx2 is actually xx, but changes in columns. u is the frequencies changing in
% %lines.

x1=x+xx/2;
x2=x-xx/2;
% 
% % equation (12) 'cross spectral density in far field'
% 
% sigma=3e-3; 
% % x1-x2 direction; inverse propotional
% % also proptionalto the u expansion in wigner
% 
% 
% theta=0.1; % x1+x2 direction; inverse propotional



% wavelength 2.25 micro, sigma=delta=5e-4

% M=40; 

k=2*pi/lambda;

% r=1e6;

C=0;  

SUM=0;

for m=1:M
%        m=M;
    
     C_m=binom_koef(M,m)*(-1)^(m-1)/sqrt(m);
%     C_m=nchoosek(M,m)*(-1)^(m-1)/sqrt(m);
    
% y=binary_coefficient(M,m);
% 
% C_m=y*(-1)^(m-1)/sqrt(m);
    
% C_m=factorial(M)/factorial(m)/factorial(M-m)*(-1)^(m-1)/sqrt(m);

C=C+C_m;

a_mx=sigma*sqrt((2*m*delta^2+4*sigma^2)/(m*delta^2+4*sigma^2));

b_mx=sqrt(2/m/delta^2+1/sigma^2);

c_mx=k^2*sigma^2*m*delta^2/(m*delta^2+4*sigma^2);

d_mx=2*k^2*sigma^4/(m*delta^2+4*sigma^2);

SUM_m=C_m*a_mx/b_mx*exp(-c_mx*(x1.^2+x2.^2)-d_mx*(x1-x2).^2);

% SUM_m=C_m*a_mx/b_mx*exp(-c_mx*(x1+x2).^2-d_mx*(x1-x2).^2);

SUM=SUM+SUM_m;

end

%

Corr_rect_GS=2*k^2*cos(theta)^2/C^2/r^2.*SUM;%.*exp(1i*k*r*(x2.^2-x1.^2));

Corr_rect_GS=Corr_rect_GS/max(max(Corr_rect_GS));


Phase=exp(-2i*pi/lambda*u.*xx2);
% corresponding-pixel multiplication 
% phase shift factor, horizontal axis is xx (rr), vertical axis is q.

Wig_rect=real(Phase*Corr_rect_GS*dx);

Wig_rect=Wig_rect/max(max(Wig_rect));


real_Corr_rect_GS=real(Corr_rect_GS);

abs_Corr_rect_GS=abs(Corr_rect_GS);

imag_Corr_rect_GS=imag(Corr_rect_GS);

plot_2d(xi,'x=(x1+x2)/2',x_prime,'x prime=x1-x2',real_Corr_rect_GS,...
    'Cross Spectral Density of Rectangular GS (real)',48,path);

% plot_2d(xi,'x=(x1+x2)/2',x_prime,'x prime=x1-x2',abs_Corr_rect_GS,...
%     'Cross Spectral Density of Rectangular GS (abs)',42,path);

% plot_2d(xi,'x=(x1+x2)/2',x_prime,'x prime=x1-x2',imag_Corr_rect_GS,...
%     'Cross Spectral Density of Rectangular GS (imag)',43,path);

% figure(42)
% surf(xi,x_prime,abs(Corr_rect_GS));shading interp;title('Cross Spectral Density of Rectangular GS (real)');

plot_2d(xi,'x (mm)',ui,'u (radian)',Wig_rect,'Wigner of Rectangular GS',46,path);



I_x=Corr_rect_GS(x_prime==0,:);

I_x_rect_GS=I_x/max(I_x);

I_u_rect_GS=sum(Wig_rect,2)/max(sum(Wig_rect,2));

plot_1d(xi,'x (mm)',I_x_rect_GS,'I','-','Rectangular GS I(x)',44,path);

plot_1d(ui,'u',I_u_rect_GS,'I(u)','-','Rectangular GS I(u)',41,path);

end

% test test
