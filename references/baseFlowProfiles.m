clear all;clc;close all

numbers=100;

%%%%%%%%%%%%%%%%%% Default [0,1] %%%%%%%%%%%%%%%%%%%%
ymin = 0;
ymax = 1;
y=[ymin:(ymax-ymin)/numbers:ymax];
f_ext=1-y;
f_apr=1-y;
figure()
plot(y,f_ext,'r','LineWidth',1.5);hold on
plot(y,f_apr,'b--','LineWidth',1.5);
xlabel('x','FontSize',13,'FontWeight','bold');
ylabel('y','FontSize',13,'FontWeight','bold');
set(gca,'fontsize',13)
title('Default Base flow','FontSize',13,'FontWeight','bold')
box
legend('Exact','Fourier')
grid on

%%%%%%%%%%%%%%%%%% Eaves2016JFM [-1,1] %%%%%%%%%%%%%%%%%%%%
ymin = -1;
ymax = 1;
y=[ymin:(ymax-ymin)/numbers:ymax];
R=20;
base_u= y;
base_rho_ext=-1/2 * (tanh(R*(y-1/3))+tanh(R*(y+1/3)));
base_rho_fourier=(R*(R^2/2 - (R^2*sinh(R/3)^2)/cosh(R/3)^2) - R^3/6 + (R^3*sinh(R/3)^2)/(2*cosh(R/3)^2) - (sinh(R/3)*((R^3*sinh(R/3))/(3*cosh(R/3)) + (R*sinh(R/3)*(R^2/2 - (R^2*sinh(R/3)^2)/cosh(R/3)^2))/cosh(R/3)))/cosh(R/3))*y.^3 + ((R*sinh(R/3)^2)/cosh(R/3)^2 - R)*y;
figure()
plot(base_u,y,'r','LineWidth',1.5);hold on
plot(base_rho_ext,y,'b','LineWidth',1.5);hold on
plot(base_rho_fourier,y,'g--','LineWidth',1.5);
xlabel('U, \rho','FontSize',13);
ylabel('y','FontSize',13);
set(gca,'fontsize',13)
title('Base flow of Eaves2016JFM','FontSize',13,'FontWeight','bold')
box
legend('U_{ext}','\rho_{ext}','\rho_{fourier}')
grid on



