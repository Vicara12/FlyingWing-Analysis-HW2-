
% -------------------------------------------------------------------------     
% PROGRAM LLWING: Weissinger lifting-line method for trapezoidal wings
% 220024 Aerodinàmica ESEIAAT-UPC
% e.ortega@upc.edu
% -------------------------------------------------------------------------     

%clc;
clear all;
close all;
format long;

% -------------------------------------------------------------------------
% INPUT DATA
% -------------------------------------------------------------------------

% Wing planform (assumes planar wing)

AR = 21.3 ;   % aspect ratio
TR = 0.18 ;   % taper ratio
DE25 = 17 ; % sweep angle at c/4 (deg)

ETIP = -7.1; % tip twist (deg, negative for washout)

% Sections data (uses linear interpolation between root and tip)

A0p = [ -0.101 -0.111 ]; % root and tip section zero-lift angles (deg)
CM0p = [ 0.003 0.003 ]; % root and tip section free moments
CDP = [ 0.00575 -0.0079 0.0137  ;   % root section CD0, k1 and k2  (airfoil CD curve)
        0.00692 -0.00767 0.014 ] ;  % tip section CD0, k1 and k2

% Flap/aileron (symmetrical deflection)

YF_pos = [ 0.0 0.0]; % 2y/b initial and final position of the flap/aileron in the half-wing
CF_ratio = 0.2 ;  % flap_chord/chord ratio
DE_flap = 0.0; % flap deflection (deg, positive:down)
FlapCorr = 1.0 ; % flap effectiviness (<=1)


% perform only convergence analysis
convergence = false;

if convergence
    N_conv = [4 8 16 32 64 128 256 518 1024 2048];
    ALPHA = [4.0];
    results = zeros(size(N_conv));

    for i = 1:size(N_conv, 2)
        N = N_conv(1,i)
        [c4nods,c75nods,chord,s_pan,h,Cm0_y,normals,mac,S] = geo(AR,TR,N,DE25,ETIP,A0p,CM0p,CDP,YF_pos,CF_ratio,DE_flap,FlapCorr); 
        [inv_A,wake_len] = infcoeff(N,c4nods,c75nods,normals,h) ;
        [GAMMA,Ui,ncases] = getcirc(N,ALPHA,inv_A,normals) ;
        [cl_local,force_coeff] = KuttaJoukowsky(N,c4nods,h,GAMMA,Ui,s_pan,Cm0_y,chord,CDP,ncases,wake_len,S,mac,ALPHA) ;
        results(1,i) = force_coeff(3, 1);
    end
    
    plot(N_conv, results); grid on;
    title('Convergence test for the number of panels');
    xlabel('number of panels');
    ylabel('force per unit chord (N/m)');
    
    return;
end

% Simulation data (by the time being only longitudinal analysis)

N = 128; % number of panels along the span

ALPHA = [ -10. -8.0 -4.0 0. 4.0 8.0 10. ]; % angles of attack for analysis (deg) 

% -------------------------------------------------------------------------
% LIFTING LINE SOLUTION
% -------------------------------------------------------------------------

% Wing discretization (lenghts are dimensionless with the wing span)

[c4nods,c75nods,chord,s_pan,h,Cm0_y,normals,mac,S] = geo(AR,TR,N,DE25,ETIP,A0p,CM0p,CDP,YF_pos,CF_ratio,DE_flap,FlapCorr); 

% Assembly of the influence coefficients matrix (needed once)

[inv_A,wake_len] = infcoeff(N,c4nods,c75nods,normals,h) ;

% Solve circulations for different angles of attack

[GAMMA,Ui,ncases] = getcirc(N,ALPHA,inv_A,normals) ;

% Loads calculation using plane Kutta-Joukowsky theorem (costlier, but general... and illustrative!)

[cl_local,force_coeff] = KuttaJoukowsky(N,c4nods,h,GAMMA,Ui,s_pan,Cm0_y,chord,CDP,ncases,wake_len,S,mac,ALPHA) ;

% -------------------------------------------------------------------------
% POSTPROCESS ( your work starts here...)
% -------------------------------------------------------------------------

c_avg = 18.8/20; % c = S/b mean aerodynamic chord
x_cg = 1.3; % position of the gravity center

%
%       SECTION 2
%   

% 2.1: Cl calculations
subplot(3, 2, 1);
Cl_slope = polyfit(ALPHA, force_coeff(7,:), 1);
alpha_l0 = -Cl_slope(2)/Cl_slope(1);
fprintf('Cl slope\n');
fprintf('  * alpha Cl0: %f dg\n  * dCl/dalpha: %f 1/dg\n\n', alpha_l0, Cl_slope(1));
hold on;
v = yline(0);
h = xline(0);
v.Color = [.4 .4 .4];
h.Color = [.4 .4 .4];
f1 = fplot(poly2sym(Cl_slope));
f2 = scatter(ALPHA, force_coeff(7,:), 20, 'filled');
grid on;
title('Cl vs alpha');
xlabel('alpha (dg)');
ylabel('Cl');
lgd = legend([f1, f2], 'regression', 'values');
lgd.Location = 'southeast';

% 2.2: Cm calculations
subplot(3, 2, 2);
Cm_slope = polyfit(force_coeff(7,:), force_coeff(5,:), 1);
X_ac = -Cm_slope(1)*c_avg;
fprintf('Cm slope\n');
fprintf('  * Cm0: %f\n  * Xac: %f \n\n', Cm_slope(2), X_ac);
hold on;
v = yline(0);
h = xline(0);
v.Color = [.4 .4 .4];
h.Color = [.4 .4 .4];
f1 = fplot(poly2sym(Cm_slope), [min(force_coeff(7,:))*1.2, max(force_coeff(7,:))*1.2]);
f2 = scatter(force_coeff(7,:), force_coeff(5,:), 20, 'filled');
grid on;
title('Cm vs Cl');
xlabel('Cl');
ylabel('Cm');
lgd = legend([f1, f2], 'regression', 'values');
lgd.Location = 'northeast';

% 2.3: basic and aditional lift calculations
i=1;
j=3;
Cla = (cl_local(:,i) - cl_local(:,j))/(force_coeff(7,i) - force_coeff(7,j));
Clb = cl_local(:,i) - Cla*force_coeff(7,i);
subplot(3,2,3);
plot(1:N, Cla, 1:N, Clb);
xlim([0, N+1]);
ylim([min(Clb)*1.2 max(Cla)*1.2]);
grid on;
title('Basic and aditional lift calculations');
xlabel(sprintf('span (from 1 to %d)', N));
ylabel('Cl');
lgd = legend('Cla', 'Clb');
lgd.Location = 'east';

% 2.4: Cd calculations
subplot(3, 2, 5);
Cd_slope = polyfit(force_coeff(7,:), force_coeff(11,:), 2);
fprintf('Cd slope\n');
fprintf('*  Cd = %f*Cl^2 + %f*Cl + %f\n\n', Cd_slope(1), Cd_slope(2), Cd_slope(3));
hold on;
v = yline(0);
h = xline(0);
v.Color = [.4 .4 .4];
h.Color = [.4 .4 .4];
f1 = fplot(poly2sym(Cd_slope), [min(force_coeff(7,:))*1.2, max(force_coeff(7,:))*1.2]);
f2 = scatter(force_coeff(7,:), force_coeff(10,:), 20, 'filled');
grid on;
title('Cd vs Cl');
xlabel('Cl');
ylabel('Cd');
lgd = legend([f1, f2], 'regression', 'values');
lgd.Location = 'northeast';

% 2.5: Cmcg calculations
subplot(3, 2, 6);
Cmcg = Cm_slope(2) - force_coeff(7,:)*(X_ac/c_avg - x_cg/0.94);
Cmcg_slope = polyfit(force_coeff(7,:), Cmcg, 1);
Cl_trim = -Cmcg_slope(2)/Cmcg_slope(1);
fprintf('Cmcg slope\n');
fprintf('*  Cmcg = %f*Cl + %f\n', Cmcg_slope(1), Cmcg_slope(2));
fprintf('*  Cl_trim: %f\n\n', Cl_trim);
hold on;
v = yline(0);
h = xline(0);
v.Color = [.4 .4 .4];
h.Color = [.4 .4 .4];
f1 = fplot(poly2sym(Cmcg_slope), [min(force_coeff(7,:))*1.2, max(force_coeff(7,:))*1.2]);
f2 = scatter(force_coeff(7,:), Cmcg, 20, 'filled');
grid on;
title('Cmcg vs Cl');
xlabel('Cl');
ylabel('Cmcg');
lgd = legend([f1, f2], 'regression', 'values');
lgd.Location = 'northeast';