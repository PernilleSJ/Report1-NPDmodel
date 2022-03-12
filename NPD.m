function dydt=NPD(t,y,param)
% Ja=zeros(1,p.n+1);
% Jd=zeros(1,p.n+1);
% J=zeros(1,p.n+1);
% dphidt=zeros(1,p.n);

P=y(1:param.n);
N=y(param.n+1:2*param.n);
D=y(2*param.n+1:end);

%% Phytoplankton
%Advenction and diffusion
%For loop to calculate interior fluxes of phytoplankton
for i = 2:param.n
    Ja_P(i)=param.v*P(i-1); %avection
    Jd_P(i)=-param.T_d*(P(i)-P(i-1))/param.dz; %diffusion
end

%Equation for boundary fluxes
Ja_P(1)=0; %Surface
Jd_P(1)=0; %Surface
Ja_P(param.n+1)=0; %Bottom
Jd_P(param.n+1)=0; %Bottom

J_P = Ja_P + Jd_P;

%Assemble
for i = 1:param.n %1 because we start at the bottom
    dPdt(i)=-((J_P(i+1)-J_P(i))/param.dz);
end


 %% Nutrients
 %Advection diffusion for nutrients:
 %For loop to calculate interior fluxes of nutrients
 for i = 2:param.n
     Ja_N(i)=0; %no advenction/sinking because nutrients don't sink
     Jd_N(i)=-param.T_d*(N(i)-N(i-1))/param.dz; %diffusion
 end

 %Equation for boundary fluxes
 Ja_N(1)=0; %Surface
 Jd_N(1)=0; %Surface
 Ja_N(param.n+1)=0; %Bottom
 Jd_N(param.n+1) = -param.T_d*((param.N_b-N(param.n))/param.dz); %Bottom
 J_N = Ja_N + Jd_N;
 
 %Assemble
 for i = 1:param.n %1 because we start at the bottom
     dNdt(i)=-((J_N(i+1)-J_N(i))/param.dz);
 end


%% Detritus
%Advenction and diffusion
%For loop to calculate interior fluxes of detritus
for i = 2:param.n
    Ja_D(i)=param.w*D(i-1); %avection
    Jd_D(i)=-param.T_d*(D(i)-D(i-1))/param.dz; %diffusion
end

%Equation for boundary fluxes
Ja_D(1)=0; %Surface
Jd_D(1)=0; %Surface
Ja_D(param.n+1)=0; %Bottom
Jd_D(param.n+1)=0; %Bottom

J_D = Ja_D + Jd_D;

%Assemble
for i = 1:param.n %1 because we start at the bottom
    dDdt(i)=-((J_D(i+1)-J_D(i))/param.dz);
end

%% Adding light 
I=lightfunc(P',param);

%% Adding Growth of phytoplankton
mu=param.mu_max*min(N./(param.H_n+N) , I./(param.H_i+I));
%mu*P is the growth rate for phytoplankton
% 
% %% The final equations 

% dPdt = dPdt+mu.*P'-param.L.*P';
% dNdt = dNdt-param.alpha*mu.*P'+param.tau.*D';
% dDdt = dDdt+(param.alpha*param.L.)*P'-param.tau.*D';
% dydt = [dPdt dNdt dDdt]';

%Maybe this instead: L*P-tau*D+dDdt

%% Assemble
%change of phyplankton
%dPdt=growth-loss-sinking+mixing
P_growth=mu.*P;
P_loss=param.L.*P;
P_sink_mix=dPdt';

dPdt=P_growth-P_loss+P_sink_mix;

%change of nutrients
%dNdt=-uptake+recycling+mixing
N_up=param.alpha.*P_growth;
N_detri=param.tau*D;
N_mix=dNdt';
dNdt=-N_up+N_mix+N_detri;

%change of detritus
D_in=P_loss;
D_out=N_detri;
D_mix=dDdt';

dDdt=D_in-D_out+D_mix;

%Result
dydt=[dPdt',dNdt',dDdt']';

end 