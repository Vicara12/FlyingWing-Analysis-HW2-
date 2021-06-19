
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

YF_pos = [ 0.0 1.0]; % 2y/b initial and final position of the flap/aileron in the half-wing
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
        results(1,i) = force_coeff(7, 1);
    end
    
    plot(log2(N_conv), results); grid on;
    title('Convergence test for the number of panels');
    xlabel('log2(number of panels)');
    ylabel('Cl');
    yline(max(results), '--', round(max(results), 5));
    ylim([min(results) max(results)*1.05]);
    xlim([1 max(log2(N_conv))*1.1]);
    
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

b = 20; % m (wing span)
s = 18.8; % m^2 (wing surface)
c_avg = s/b; % c = S/b mean aerodynamic chord
x_cg = 1.3; % position of the gravity center
ws = 20; % kg/m^2 wing load (W/s)
rho = 1.225; % kg/m^3 air density
v_flight = 25; % m/s aircraft's design speed
Clmax_root = 1.3122; % stall Cl for the wing root (at 12 dg)
Clmax_tip = 1.2465; % stall Cl for the wing tip (at 12 dg)


%
%       SECTION 2
%   

% 2.1: Cl calculations
subplot(3, 3, 1);
Cl_slope = polyfit(ALPHA, force_coeff(7,:), 1);
alpha_l0 = -Cl_slope(2)/Cl_slope(1);
fprintf('\n\n2.1: Cl slope\n');
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
subplot(3, 3, 2);
Cm_slope = polyfit(force_coeff(7,:), force_coeff(5,:), 1);
X_ac = -Cm_slope(1)*c_avg;
fprintf('2.2: Cm slope\n');
fprintf('  * Cm = %f*Cl + %f\n', Cm_slope(1), Cm_slope(2));
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
subplot(3,3,3);
plot(1:N, Cla, 1:N, Clb);
xlim([0, N+1]);
ylim([min(Clb)*1.2 max(Cla)*1.2]);
grid on;
title('2.3: Basic and aditional lift calculations');
xlabel(sprintf('span (from 1 to %d)', N));
ylabel('Cl');
lgd = legend('Cla', 'Clb');
lgd.Location = 'east';


% 2.4: Cd calculations
subplot(3, 3, 4);
Cd_slope = polyfit(force_coeff(7,:), force_coeff(11,:), 2);
fprintf('2.4: Cd slope\n');
fprintf('  * Cd = %f*Cl^2 + %f*Cl + %f\n\n', Cd_slope(1), Cd_slope(2), Cd_slope(3));
hold on;
v = yline(0);
h = xline(0);
v.Color = [.4 .4 .4];
h.Color = [.4 .4 .4];
f1 = fplot(poly2sym(Cd_slope), [min(force_coeff(7,:))*1.2, max(force_coeff(7,:))*1.2]);
f2 = scatter(force_coeff(7,:), force_coeff(11,:), 20, 'filled');
grid on;
title('Cd vs Cl');
xlabel('Cl');
ylabel('Cd');
lgd = legend([f1, f2], 'regression', 'values');
lgd.Location = 'northeast';


% 2.5: Cmcg calculations
subplot(3, 3, 5);
Cmcg = Cm_slope(2) - force_coeff(7,:)*(X_ac/c_avg - x_cg/0.94);
Cmcg_slope = polyfit(force_coeff(7,:), Cmcg, 1);
Cl_trim = -Cmcg_slope(2)/Cmcg_slope(1);
fprintf('2.5: Cmcg slope\n');
fprintf('  * Cmcg = %f*Cl + %f\n', Cmcg_slope(1), Cmcg_slope(2));
fprintf('  * Cl_trim: %f\n\n', Cl_trim);
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




%
%       SECTION 3
%


% 3.1: search for wing sweep that gives 15% distance
Cm_slope = polyfit(force_coeff(7,:), force_coeff(5,:), 1);
X_ac = -Cm_slope(1)*c_avg;

fprintf('\n3.1: Wing sweep\n');
fprintf('  * Current distance cg ac: %.2f%%\n', 100*(X_ac/c_avg - x_cg/c_avg));
margin = 0.001;
step = 1;
desired_Xac = c_avg*(x_cg/c_avg + .15);
last_positive = true;
fprintf('\nComputing optimal wing sweep...\n\n');

while (abs(desired_Xac - X_ac) > margin)
    
    [c4nods,c75nods,chord,s_pan,h,Cm0_y,normals,mac,S] = geo(AR,TR,N,DE25,ETIP,A0p,CM0p,CDP,YF_pos,CF_ratio,DE_flap,FlapCorr); 
    [inv_A,wake_len] = infcoeff(N,c4nods,c75nods,normals,h) ;
    [GAMMA,Ui,ncases] = getcirc(N,ALPHA,inv_A,normals) ;
    [cl_local,force_coeff] = KuttaJoukowsky(N,c4nods,h,GAMMA,Ui,s_pan,Cm0_y,chord,CDP,ncases,wake_len,S,mac,ALPHA) ;
    
    Cm_slope = polyfit(force_coeff(7,:), force_coeff(5,:), 1);
    X_ac = -Cm_slope(1)*c_avg;
    
    difference = desired_Xac - X_ac;
    
    if (difference > 0)
        if (~last_positive)
            step = step/2;
        end
        DE25 = DE25 + step;
    else
        if (last_positive)
            step = step/2;
        end
        DE25 = DE25 - step;
    end
                
    last_positive = difference > 0;
end

fprintf('  * New sweep angle: %.2fº\n  * New Xac: %f\n  * Distance ac cg: %.2f%%\n', DE25, X_ac, 100*(X_ac - x_cg)/c_avg);


% 3.2: new wing's geometric torsion
Cl_des = 2*9.81*ws/(rho*v_flight^2);
Cmcg = Cm_slope(2) - force_coeff(7,:)*(X_ac/c_avg - x_cg/0.94);
Cmcg_slope = polyfit(force_coeff(7,:), Cmcg, 1);
Cl_trim = -Cmcg_slope(2)/Cmcg_slope(1);
fprintf('\n3.2: Wing washout\n');
fprintf('  * Current trimmed Cl: %f\n  * Desired trimmed Cl: %f\n', Cl_trim, Cl_des);
margin = 0.0001;
step = 1;
last_positive = true;
fprintf('\nComputing optimal wing washout...\n\n');

while (abs(Cl_des - Cl_trim) > margin)
    
    [c4nods,c75nods,chord,s_pan,h,Cm0_y,normals,mac,S] = geo(AR,TR,N,DE25,ETIP,A0p,CM0p,CDP,YF_pos,CF_ratio,DE_flap,FlapCorr); 
    [inv_A,wake_len] = infcoeff(N,c4nods,c75nods,normals,h) ;
    [GAMMA,Ui,ncases] = getcirc(N,ALPHA,inv_A,normals) ;
    [cl_local,force_coeff] = KuttaJoukowsky(N,c4nods,h,GAMMA,Ui,s_pan,Cm0_y,chord,CDP,ncases,wake_len,S,mac,ALPHA) ;
    
    Cm_slope = polyfit(force_coeff(7,:), force_coeff(5,:), 1);
    Cmcg = Cm_slope(2) - force_coeff(7,:)*(X_ac/c_avg - x_cg/0.94);
    Cmcg_slope = polyfit(force_coeff(7,:), Cmcg, 1);
    Cl_trim = -Cmcg_slope(2)/Cmcg_slope(1);
    
    difference = Cl_des - Cl_trim;
    
    %fprintf('difference: %f      Cl_des: %f       Cl_trim: %f       washout: %f\n', difference, Cl_des, Cl_trim, ETIP);
    
    if (difference > 0)
        if (~last_positive)
            step = step/2;
        end
        ETIP = ETIP - step;
    else
        if (last_positive)
            step = step/2;
        end
        ETIP = ETIP + step;
    end
                
    last_positive = difference > 0;
end

fprintf('  * New washout angle: %.2fº\n  * New Cl_trim: %f\n', ETIP, Cl_trim);


% 3.3: Stall characteristics of the wing
DE25 = 17.75;
ETIP = -4.66;
[c4nods,c75nods,chord,s_pan,h,Cm0_y,normals,mac,S] = geo(AR,TR,N,DE25,ETIP,A0p,CM0p,CDP,YF_pos,CF_ratio,DE_flap,FlapCorr); 
[inv_A,wake_len] = infcoeff(N,c4nods,c75nods,normals,h) ;
[GAMMA,Ui,ncases] = getcirc(N,ALPHA,inv_A,normals) ;
[cl_local,force_coeff] = KuttaJoukowsky(N,c4nods,h,GAMMA,Ui,s_pan,Cm0_y,chord,CDP,ncases,wake_len,S,mac,ALPHA) ;

Clmax_slope = polyfit([0 b/2], [Clmax_root Clmax_tip], 1);
Clmax_values = polyval(Clmax_slope, (1:(N/2))/N*b);
Clmax_values = [fliplr(Clmax_values) Clmax_values]';
i=1;
j=3;
Cla = (cl_local(:,i) - cl_local(:,j))/(force_coeff(7,i) - force_coeff(7,j));
Clb = cl_local(:,i) - Cla*force_coeff(7,i);
CL = (Clmax_values - Clb)./Cla;
span_values = (1:(N/2))/N*b;
span_values = [-fliplr(span_values) span_values]';
subplot(3, 3, 6)
hold on;
p1 = plot(span_values, Clmax_values);
p2 = plot(span_values, CL);
grid on;
xline(b/2*0.7, '--', 'control surfaces');
xline(-b/2*0.7, '--', 'control surfaces');
title('CL and Clmax along the uncorrected wing span');
xlabel('wing span (m)');
ylabel('Cl');
lgd = legend([p1 p2], 'Cl max', 'CL');
lgd.Location = 'north';
%{
DE25 = 17.75;
ETIP = -2.66;
[c4nods,c75nods,chord,s_pan,h,Cm0_y,normals,mac,S] = geo(AR,TR,N,DE25,ETIP,A0p,CM0p,CDP,YF_pos,CF_ratio,DE_flap,FlapCorr); 
[inv_A,wake_len] = infcoeff(N,c4nods,c75nods,normals,h) ;
[GAMMA,Ui,ncases] = getcirc(N,ALPHA,inv_A,normals) ;
[cl_local,force_coeff] = KuttaJoukowsky(N,c4nods,h,GAMMA,Ui,s_pan,Cm0_y,chord,CDP,ncases,wake_len,S,mac,ALPHA) ;

i=1;
j=3;
Cla = (cl_local(:,i) - cl_local(:,j))/(force_coeff(7,i) - force_coeff(7,j));
Clb = cl_local(:,i) - Cla*force_coeff(7,i);
CL = (Clmax_values - Clb)./Cla;
subplot(3, 3, 7)
hold on;
p1 = plot(span_values, Clmax_values);
p2 = plot(span_values, CL);
grid on;
xline(b/2*0.7, '--', 'control surfaces');
xline(-b/2*0.7, '--', 'control surfaces');
title('CL and Clmax along the corrected wing span');
xlabel('wing span (m)');
ylabel('Cl');
lgd = legend([p1 p2], 'Cl max', 'CL');
lgd.Location = 'north';
%}
v_stall = sqrt(2*ws*9.81/(min(CL)*rho));
fprintf('\n3.3: Wing stall characteristics\n');
fprintf('  * Stall speed: %.2f m/s\n', v_stall);


% 3.4: wing's efficiency
DE25 = 17.75;
ETIP = -4.66;
[c4nods,c75nods,chord,s_pan,h,Cm0_y,normals,mac,S] = geo(AR,TR,N,DE25,ETIP,A0p,CM0p,CDP,YF_pos,CF_ratio,DE_flap,FlapCorr); 
[inv_A,wake_len] = infcoeff(N,c4nods,c75nods,normals,h) ;
[GAMMA,Ui,ncases] = getcirc(N,ALPHA,inv_A,normals) ;
[cl_local,force_coeff] = KuttaJoukowsky(N,c4nods,h,GAMMA,Ui,s_pan,Cm0_y,chord,CDP,ncases,wake_len,S,mac,ALPHA);
Cl_slope = polyfit(ALPHA, force_coeff(7,:), 1);
Cd_slope = polyfit(force_coeff(7,:), force_coeff(11,:), 2);
speed = round(v_stall):50;
cl = 2*ws*9.81./(rho*speed.*speed);
E = 1./(Cd_slope(3)./cl + Cd_slope(2) + Cd_slope(1).*cl);
Cd0 = Cd_slope(3);
k2 = Cd_slope(1);
subplot(3, 3, 8);
hold on;
plot(speed, E);
title('Aerodynamic efficiency vs speed');
xlabel('speed (m/s)');
ylabel('efficiency (Cl/Cd)');
grid on;
grid minor;
xlim([v_stall-3 50+3]);
fprintf('\n3.4: Wing efficiency analysis\n');
fprintf('  * 1/E = %f/Cl + %f + %f*Cl\n', Cd_slope(3), Cd_slope(2), Cd_slope(1));
fprintf('  * Max efficiency: %f\n  * Speed for max efficiency: %f\n', max(E), speed(E == max(E)));
fprintf('  * Cl for maximum efficiency: %f\n', cl(E == max(E)));


% 3.5: Flap deflection analysis
flap_angles = linspace(-20, 20, 60);
speed = zeros(size(flap_angles));
cl_trim = zeros(size(flap_angles));
fprintf('\n3.5: Flap deflection analysis\n');

for i = 1:size(flap_angles, 2)
    
    DE_flap = flap_angles(i);
    
    [c4nods,c75nods,chord,s_pan,h,Cm0_y,normals,mac,S] = geo(AR,TR,N,DE25,ETIP,A0p,CM0p,CDP,YF_pos,CF_ratio,DE_flap,FlapCorr); 
    [inv_A,wake_len] = infcoeff(N,c4nods,c75nods,normals,h) ;
    [GAMMA,Ui,ncases] = getcirc(N,ALPHA,inv_A,normals) ;
    [cl_local,force_coeff] = KuttaJoukowsky(N,c4nods,h,GAMMA,Ui,s_pan,Cm0_y,chord,CDP,ncases,wake_len,S,mac,ALPHA) ;

    Cm_slope = polyfit(force_coeff(7,:), force_coeff(5,:), 1);
    X_ac = -Cm_slope(1)*c_avg;
    Cmcg = Cm_slope(2) - force_coeff(7,:)*(X_ac/c_avg - x_cg/0.94);
    Cmcg_slope = polyfit(force_coeff(7,:), Cmcg, 1);
    cl_trim(i) = -Cmcg_slope(2)/Cmcg_slope(1);
    if (cl_trim(i) > 0)
        speed(i) = sqrt(2*ws*9.81/(rho*cl_trim(i)));
    end
end

subplot(3,3,9);
hold on;
yyaxis left;
p1 = plot(flap_angles(cl_trim > 0 & speed < 100), speed(cl_trim > 0 & speed < 100));
grid on;
title('Trimmed leveled flight conditions vs flap deflection');
xlabel('flap deflection angle (dg)');
ylabel('trimmed leveled flight speed (m/s)');
yyaxis right;
p2 = plot(flap_angles, cl_trim);
ylabel('trimmed leveled flight Cl');

Cl_emax = sqrt(Cd0/k2);
Cl_mindes = sqrt(Cd0/3*k2);
flap_angle_to_cl = polyfit(cl_trim, flap_angles, 10);
flap_cl_emax = polyval(flap_angle_to_cl, Cl_emax);
flap_cl_mindes = polyval(flap_angle_to_cl, Cl_mindes);

scatter([0.142916 6.506799],[0.487613 0.008066], 30, 'filled', 'g');
lgd = legend([p1 p2], 'spped', 'Cl');

fprintf('  * Cl max efficiency: %f\n  * Flap deflection for Cl max e: %f dg\n', Cl_emax, flap_cl_emax);
fprintf('  * Cl min discent rate: %f\n  * Flap deflection for Cl min descent rate: %f dg\n', Cl_mindes, flap_cl_mindes);

fprintf('\n\nPROGRAM HAS FINISHED\n\n');