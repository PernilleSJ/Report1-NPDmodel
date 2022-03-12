% ======================================================
% Run a advection-diffusion equation
% ======================================================
function [t,z,P] = grid_func(dz) % dz is the grid spacing
%
% Parameters:
%
param.T_d = 4; % Diffusivity (m2/day)
param.v = 0.96; % Settling velocity (m/day)

param.depth = 100; % meter
param.dz = dz; % Grid spacing (m)
param.z = 0.5*param.dz:param.dz:(param.depth-0.5*param.dz);
param.nGrid = length(param.z);  % No. of grid cells
%
% Initialization:
%
P0 = exp(-(param.z-100).^2/5);
%
% Run model
%
[t, P] = ode45(@PmodelDeriv, [0:200], P0, [], param);
z = param.z;


    % ---------------------------------------------------------
    % Derivative function
    % ---------------------------------------------------------
    function dPdt = PmodelDeriv(t,P,param)
        %
        % Advective fluxes
        %
        i = 2:(param.nGrid);
        Ja_P(i) = param.v*P(i-1);
        Ja_P(1) = 0;  % No input from the surface
        Ja_P(param.nGrid+1) = 0; % Closed bottom
        %
        % Diffusive fluxes:
        %
        Jd_P(i) = -param.T_d*(P(i)-P(i-1))/param.dz;
        Jd_P(1) = 0; % No flux at the surface...
        Jd_P(param.nGrid+1) = 0;  % ...or the bottom
        %
        % Rate-of-change due to advection and diffusion:
        %
        J = Ja_P + Jd_P;
        dPdt = -(J(2:(param.nGrid+1))-J(1:param.nGrid))/param.dz;
        % Make dPdt a column vector:
        dPdt = dPdt';
    end

end