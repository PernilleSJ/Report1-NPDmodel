%% NPD model
clear all
close all
clc

%% -------Set parameters-------
param.v= 1.008;           %sinking velocity [m/day]
param.v0= 0;              %sinking velocity [m/day]
param.v2= 2;              %sinking velocity [m/day]
param.T_d= 10;            %turbulent diffusion rate [m^2/day] 
param.depth= 100;         %depth of water column [m]
param.n= 100;             %number of grid cells
param.I0=350;             % incident light intensity [mymol photons m^-2 s^-1]
param.Kbg=0.045;          %background turbidity [1/m]
param.k=15*10^-12;        %specific light attenuation of phytoplankton [m^2/cell]
param.H_i=30;             %Half satuaration ligth growth [mumol photons/(m^2/s)]
param.pmax=0.96;          %%maximum specific production rate  [1/day]
param.L=0.24;             %specific loss rate [1/day]
param.alpha=1*10^-9;      %nutrient content of pp [mmol N/(m^2s)]
param.eps=0.03;           % (1/day)
param.gamma=0.01;         % (1/day)
param.N_b=50;             %bottom concentration af nutrients [mmol/m^3]
param.mu_max=0.96;        %max speciffic growth rate [1/d]
param.H_n=0.02;  %0.0425;         %Half satuaration nutrients growth [mmol N/(m^3)]
%Detritus parameters
param.tau=0.1;            %Remineralizaytion 1/d
param.w=5;                %sinking speed of detritus m/d


%% Making grid
param.dz=param.depth/param.n; %width of grid
param.z=0.5*param.dz:param.dz:(param.depth-0.5*param.dz); %The grid
z=0.5*param.dz:param.dz:(param.depth-0.5*param.dz); %To call it in plots

%% Initial conditions 
phi=zeros(1,param.n);
nut=zeros(1,param.n);

phif=exp(-(param.z-100).^2/5);

P0 = 2e6*exp(-(param.z-param.depth/4).^2/1000); %Gauss distribution
N0 = param.N_b*exp(-(param.z-param.depth/1.8).^2/500); %Gauss distribution
D0=zeros(1,param.n);

%Initial conditions together
y=[P0,N0,D0];


%% --------Call ode45---------
option = odeset('NonNegative',1:3*param.n);
[t,y]=ode45(@NPD, [0,2500], y, option, param);

%Define for plot
P=y(:,1:param.n);
N=y(:,param.n+1:param.n*2);
D=y(:,param.n*2+1:end);

%Light
I=lightfunc(P(end,:),param);
I_0=lightfunc(0*P(end,:),param);

%% Plots
%To see when converged
% figure()
% plot(t,N)
% xlabel('Time [days]')
% ylabel('Concentration [mmol/m^3]')
% title('Nutrient concentration in time')

%%
% figure()
% subplot(1,3,1)
% % 
% plot(I,-z,'g','Linewidth',2)
% hold on
% plot(I_0,-z,'y','Linewidth',2)
% ylabel("Depth [m]")
% xlabel("I[mumol/(m^2s^1)]")
% hold off
% title("Ligth intensity")
% legend("With plankton","Without plankton")
% grid on
% % 
% subplot(1,3,2)
% % 
% plot(N(end,:),-z,'r','Linewidth',2)
% ylabel("Depth [m]")
% xlabel("Concentration [mmol/m^3]")
% title("End distribution of Nutrients")
% grid on
% % 
% % 
% subplot(1,3,3)
% % 
% plot(P(end,:),-z,'g','Linewidth',2)
% ylabel("Depth [m]")
% xlabel("Concentration [cell/m^3]")
% title("End distribution of Phytoplankton")
% grid on
%%
%Plankton
% figure()
% subplot(2,2,1)
% % 
% surface(t,param.z,P')
% shading interp
% xlabel('time [days]')
% ylabel('depth [m]')
% title('Phytoplankton in water column')
% colorbar
% set(gca, 'YDir','reverse')
% hcb=colorbar;
% hcb.Label.String='Phytoplankton(ug m^-^3)';
% grid on
% % 
% %Nutrients
% subplot(2,2,2)
% % 
% surface(t,param.z,N')
% shading interp
% xlabel('time [days]')
% ylabel('depth [m]')
% title('Nutrients in water column')
% colorbar
% set(gca, 'YDir','reverse')
% hcb=colorbar;
% hcb.Label.String='Nutrients [mmol N/m^3]';
% grid on
% % 
% %Detritus
% % 
% subplot(2,2,3)
% surface(t,param.z,D')
% shading interp
% xlabel('time [days]')
% ylabel('depth [m]')
% title('Detritus in water column')
% colorbar
% set(gca, 'YDir','reverse')
% hcb=colorbar;
% hcb.Label.String='Detritus [mmolN/m^3]';
% grid on



 %%
% figure()
% plot(I*(1*10^9),-param.z,'y','Linewidth',2)
% hold on
% plot((1*10^-1)*N(end,:),-param.z,'r','Linewidth',2)
% hold on
% plot(P(end,:),-param.z,'g','LineWidth',2)
% ylabel("Depth [m]")
% xlabel("Concentration of I [mumol/(m^2s^1)], N [mmolN/m^3] and P [ug m^-^3]")
% hold off
% title("Light, nutrients and phytoplankton")
% legend("Light 10e9","Nutrients 10e-1",'Phytoplankton')
% xlim([0 0.7*10^11]);

%%
%Sensitivity
%subplot(1,3,3)
plot(P(end,:),-z,'b','Linewidth',2)
ylabel("Depth [m]")
xlabel("Concentration [cell/m^3]")
title('Phytoplankton with v-value = 2')

