% This code numerically solves the global sea ice and climate model
% described by Eqs. (2), (8) and (9) in Wagner & Eisenman (2015, hereafter
% WE15; see reference below). For computational efficiency, a diffusive
% 'ghost layer' is invoked, as described in WE15 Appendix A. This allows us
% to solve the system using an Implicit Euler timestep in the ghost layer
% (where diffusion occurs) and a Forward Euler timestep for the evolution
% of surface enthalpy. The document WE15_NumericIntegration.pdf provides
% further detailed documentation to go with this script.
%
% This code uses the default parameter values from WE15 (Table 1 and
% Section 2d) with no climate forcing (F=0), and the final part of the code 
% produces a plot similar to WE15 Figure 2.
%
% By default this code runs a simulation for 30 years with 1000 timesteps/year using
% a spatial resolution of 100 gridboxes, equally spaced in x=sin(lat), between the equator
% and pole.
%
% Till Wagner (tjwagner@ucsd.edu) & Ian Eisenman (eisenman@ucsd.edu), 
% created Mar 2015, minor edits Apr 2016, minor bug fix Jan 2022 [in Eq. (A1), S(:,i) -> S(:,i+1)].
%
% Reference: T.J.W. Wagner and I. Eisenman (2015). How climate model
% complexity influences sea ice stability. J Climate 28, 3998-4014.
%
%--------------------------------------------------------------------------
function [tfin_out, Efin_out] = sea_ice_EBM_WE15
%%Model parameters (WE15, Table 1 and Section 2d) -------------------------
D  = 0.6;     %diffusivity for heat transport (W m^-2 K^-1)
A  = 193;     %OLR when T = T_m (W m^-2)
B  = 2.1;     %OLR temperature dependence (W m^-2 K^-1)
cw = 9.8;     %ocean mixed layer heat capacity (W yr m^-2 K^-1)
S0 = 420;     %insolation at equator  (W m^-2)
S1 = 338;     %insolation seasonal dependence (W m^-2)
% S1 = 0;
S2 = 240;     %insolation spatial dependence (W m^-2)
a0 = 0.7;     %ice-free co-albedo at equator
a2 = 0.1;     %ice-free co-albedo spatial dependence
ai = 0.4;     %co-albedo where there is sea ice
Fb = 4;       %heat flux from ocean below (W m^-2)
k  = 2;       %sea ice thermal conductivity (W m^-2 K^-1)
Lf = 9.5;     %sea ice latent heat of fusion (W yr m^-3)
F  = 0;       %radiative forcing (W m^-2)
cg = 0.01*cw; %ghost layer heat capacity(W yr m^-2 K^-1)
tau = 1e-5;   %ghost layer coupling timescale (yr)
%%The default run in WE15 Fig 2 uses these time-stepping parameters: ------
% n  = 400;   % # of evenly spaced latitudinal gridboxes (equator to pole)
% nt = 1e3;   % # of timesteps per year (limited by numerical stability)
% dur= 200;   % # of years for the whole run
%%For a quicker computation, here we use these parameters: ----------------
% n  = 100;   % # of evenly spaced latitudinal gridboxes (equator to pole)
% nt = 1e3;   % # of timesteps per year (limited by numerical stability)
% dur= 30;    % # of years for the whole run
%%For Alignment with my results
n = 100;
nt = 2e3;
dur = 30;
dt = 1/nt;
%%Spatial Grid ------------------------------------------------------------
dx = 1/n;               %grid box width
x = (dx/2:dx:1-dx/2)';  %grid
%%Diffusion Operator (WE15 Appendix A) -----------------------------------
xb = (dx:dx:1.0-dx)';  
lambda=D/dx^2*(1-xb.^2); L1=[0; -lambda]; L2=[-lambda; 0]; L3=-L1-L2;
diffop = - diag(L3) - diag(L2(1:n-1),1) - diag(L1(2:n),-1);
%%Definitions for implicit scheme for Tg
cg_tau = cg/tau;
dt_tau = dt/tau;
dc = dt_tau*cg_tau;
kappa = (1+dt_tau)*eye(n)-dt*diffop/cg;
%%Seasonal forcing [WE15 Eq. (3)]
ty = dt/2:dt:1-dt/2;
S=repmat(S0-S2*x.^2,[1,nt])-repmat(S1*cos(2*pi*ty),[n,1]).*repmat(x,[1,nt]);
S=[S S(:,1)];
%%Further definitions
M = B+cg_tau;
aw= a0-a2*x.^2;   %open water albedo
kLf = k*Lf;
%%Initial conditions ------------------------------------------------------
T = 7.5+20*(1-2*x.^2);
% T = -x;
Tg = T; E = cw*T;
%%Set up output arrays, saving 100 timesteps/year
E100 = zeros(n,dur*100); T100 = E100;
%%Integration (see WE15_NumericIntegration.pdf)----------------------------
% Loop over Years ---------------------------------------------------------
for years = 1:dur
    % Loop within One Year-------------------------------------------------
    for i = 1:nt
        if mod(i,nt/100)==0 %store 100 timesteps per year
            E100(:,(years-1)*100+i/(nt/100)) = E;
            T100(:,(years-1)*100+i/(nt/100)) = T;
        end
        % forcing
        alpha = aw.*(E>0) + ai*(E<0);        % WE15 Eq. (4)
        C = alpha.*S(:,i) + cg_tau*Tg-A+F;  
        % surface temperature
        T0 =  C./(M-kLf./E);                 %WE15 Eq. (A3)
        T = E/cw.*(E>=0)+T0.*(E<0).*(T0<0);  %WE15 Eq. (9)
        % Forward Euler for E
        E = E+dt*(C-M*T+Fb);                 %WE15 Eq. (A2)
        % Implicit Euler for Tg
        Tg = (kappa-diag(dc./(M-kLf./E).*(T0<0).*(E<0)))\ ...
            (Tg + (dt_tau*(E/cw.*(E>=0)+(ai*S(:,i+1) ...
            -A+F)./(M-kLf./E).*(T0<0).*(E<0))));        %WE15 Eq. (A1)
    end
    if mod(years,10)==0, disp(['year ' num2str(years) ' complete']), end
end
% -------------------------------------------------------------------------
%%output only final year 
tfin = linspace(0,1,100);  
Efin = E100(:,end-99:end); 
Tfin = T100(:,end-99:end);
if nargout>0 %save output
    tfin_out=tfin;
    Efin_out=Efin;
end
% -------------------------------------------------------------------------
%WE15 Figure 2: Default Steady State Climatology --------------------------
% -------------------------------------------------------------------------
winter=26;    %time of typically coldest hemisphere-mean T
summer=76;    %time of typically warmest hemisphere-mean T
%%compute seasonal ice edge
xi = zeros(1,100);
for j = 1:length(tfin)
    Ej = Efin(:,j);
    xe=x(find(Ej<0,1,'first'));
    if ~isempty(xe)
        xi(j) = xe;
    else
        xi(j) = max(x);
    end
end
fig = figure(2); clf
set(fig, 'Position', [2561 361 1920 984]);
%%plot the enthalpy (Fig 2a)
subplot(1,4,1)
clevs = [-40:20:0 50:50:300];
contourf(tfin,x,Efin,clevs)
%%plot ice edge on E
hold on
plot(tfin,xi,'k')
colorbar
%%alternatively, use 
% cbarf(Efin,clevs);  %makes nice discrete colorbar
%%(http://www.mathworks.com/matlabcentral/fileexchange/14290-cbarf)
xlabel('t (final year)')
ylabel('x')
title('E (J/m^2)')
% plot the temperature (Fig 2b)
clevs = -30:5:30;
subplot(1,4,2)
contourf(tfin,x,Tfin,clevs)
% plot ice edge on T
hold on
plot(tfin,xi,'k')
%%plot T=0 contour (ice surface melt occurs in the region between ice edge
%%and T=0 contour)
contour(tfin,x,Tfin,[0,0],'r')
colorbar
% cbarf(Tfin,clevs);
xlabel('t (final year)')
ylabel('x')
title('T (^oC)')
%%plot the ice thickness (Fig 2c)
hfin = -(Efin/Lf.*(Efin<0));
subplot(1,4,3)
clevs = 0.0001:.5:4;
contourf(tfin,x,hfin,clevs)
%%plot ice edge on h
hold on
contour(tfin,x,hfin,[0,0])
plot([tfin(winter) tfin(winter)], [0 max(x)],'k')
plot([tfin(summer) tfin(summer)], [0 max(x)],'k--')
xlabel('t (final year)')
ylabel('x')
title('h (m)')
colorbar
% cbarf(hfin,round(clevs,1));
%%plot temperature profiles (Fig 2d)
subplot(4,4,4)
plot(x,Tfin(:,summer),'k--')
hold on
plot(x,Tfin(:,winter),'k')
plot([0 1], [0 0],'k')
xlabel('x')
ylabel('T (^oC)')
legend('summer','winter','location','southwest')
%%plot ice thickness profiles (Fig 2e)
subplot(4,4,8)
plot(x,hfin(:,summer),'k--')
hold on
plot(x,hfin(:,winter),'k')
plot([0 1], [0 0],'k')
xlim([0.7 1])
xlabel('x')
ylabel('h (m)')
%%plot seasonal thickness cycle at pole (Fig 2f)
subplot(4,4,12)
plot(tfin,hfin(end,:),'k')
xlabel('t (final year)')
ylabel('h_p (m)')
ylim([2 1.1*max(hfin(end,:))])
%%plot ice edge seasonal cycle (Fig 2g)
subplot(4,4,16)
xideg = rad2deg(asin(xi));
plot(tfin,xideg,'k-')
ylim([0 90])
xlabel('t (final year)')
ylabel('\theta_i (deg)')
% exportgraphics(gcf, "WE15.pdf", ContentType="vector")

