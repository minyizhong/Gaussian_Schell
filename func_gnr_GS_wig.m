function [Corr,Wig_GS]=func_gnr_GS_wig(nx,nu,max_x,max_u,wo,z,lambda,lc,path)


% nx is the amount of lines/columns. maximum value of x
xi=gnr_array(nx,max_x);

ui=gnr_array(nu,max_u);

% xi( xi<=-1.106 | xi>=1.094 )=0;

x_prime=xi*2;

[x,xx] = meshgrid(xi, x_prime);%(xi,xxi);
% % x changes in columns, xx changes in lines
[xx2,u]=meshgrid(x_prime,ui);%(xxi,ui);
% %xx2 is actually xx, but changes in columns. u is the frequencies changing in
% %lines.

x1=x+xx/2;
x2=x-xx/2;

dx=abs(xi(1)-xi(2));

% ------------------------------------------------------------------------------------------------------
% wo=9; % reference beam waist 
% (canceled) A=0.5; % normalized amplitude constant 
% z=100; % lateral distance from the beam waist
% lambda=500e-6; % wavelength in millimeters
% lc=7; % lateral coherent length in millimeters
k=2*pi/lambda;

wc=sqrt((lc*wo).^2/(lc.^2+wo.^2));


w_z=wo*sqrt(1+(lambda*z/pi/wo/wc).^2);
% radius of curvature at the wavefront at distance z:

if z~=0
    R_z=z*(1+(pi*wo*wc/lambda/z).^2);
    Corr=(wo/w_z).^2*exp(-(x1+x2).^2/2/w_z.^2).*exp(-wo.^2*(x1-x2).^2/2/lc.^2/w_z.^2).*exp(-1i*k*(x1.^2-x2.^2)/2/R_z);
    
else % if z=0, radius R_z is 0, at the beam waist position
%     R_z=0; 
    Corr=(wo/w_z).^2*exp(-(x1+x2).^2/2/w_z.^2).*exp(-wo^2*(x1-x2).^2/2/lc.^2/w_z.^2);
end

Corr=Corr/max(max(Corr));

Phase=exp(-2i*pi/lambda*u.*xx2);

Wig_GS=real(Phase*Corr*dx);

% --------------------------------------------------------------------------
% plot all the figures

real_Corr=real(Corr);

plot_2d(xi,'x=(x_1+x_2)/2 (mm)',x_prime,'x prime=x_2-x_1 (mm)',real_Corr,'GS correlation (real values)',2,path);


I_x_GS=sum(Wig_GS)/max(sum(Wig_GS));
plot_1d(xi,'x (mm)',I_x_GS,'I(x)','-','GS I(x)',3,path);

I_u_GS=sum(Wig_GS,2)/max(sum(Wig_GS,2));
plot_1d(ui,'u (radian)',I_u_GS,'I(u)','-', 'GS I(u)',4,path);
% % 
% plot_2d(xi,'x=(x_1+x_2)/2 (mm)',xi*2,'x prime=x_2-x_1 (mm)',imag(Corr),'GS correlation (imag values)',5,path);
% % 
% plot_2d(xi,'x=(x_1+x_2)/2 (mm)',xi*2,'x prime=x_2-x_1 (mm)',abs(Corr),'GS correlation (abs values)',6,path);
Wig_GS=Wig_GS./max(max(Wig_GS));
plot_2d(xi,'x (mm)',ui,'u (radian)',Wig_GS,'Wigner of GS beam',7,path);

% figure(1)
% pcolor(xi,ui,Wig_GS);
% colormap jet, shading flat
% colorbar('location','eastoutside')
% set(gca,'FontSize',14,'FontSize',20);
% title('Wigner of GS beam','FontSize',20)
% xlabel('x (mm)','FontSize',20)
% ylabel('u (radian)','FontSize',20)


end
