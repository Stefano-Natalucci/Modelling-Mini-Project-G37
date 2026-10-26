clear
close all
clc

%% =========================================================
%  GEOMETRY
% ==========================================================

Lbed = 0.40;

xc = Lbed/2;
yc = Lbed/2;

a = 0.15;
b = 0.10;

omega = 2*pi/10;      % one ellipse in 10 s

%% =========================================================
%  MECHANICAL PARAMETERS
% ==========================================================

m     = 0.15;
Jm    = 1.5e-6;
rm    = 6e-3;
Km    = 0.025;

mu_r  = 5.5e-3;
mu_lx = 0.04;
mu_ly = 0.04;

Lm = 3e-4;
Rm = 2;

Meq = m + 2*Jm/rm^2;

Bx = 2*mu_lx + 8*mu_r;
By = 2*mu_ly + 4*mu_r;

%% =========================================================
%  PRINT-HEAD STATE-SPACE MODEL
% ==========================================================

A = [0  1                  0                  0  0                  0;
     0 -Bx/Meq             Km/(rm*Meq)        0  0                  0;
     0 -2*Km/(Lm*rm)      -Rm/Lm              0  0                  0;
     0  0                  0                  0  1                  0;
     0  0                  0                  0 -By/Meq             Km/(rm*Meq);
     0  0                  0                  0 -2*Km/(Lm*rm)      -Rm/Lm];

B = [0       0;
     0       0;
     1/Lm    0;
     0       0;
     0       0;
     0       1/Lm];

Cmat = [1 0 0 0 0 0;
        0 0 0 1 0 0];

D = zeros(2,2);

sys = ss(A,B,Cmat,D);

%% =========================================================
%  NOZZLE PARAMETERS
% ==========================================================

Te = 25 + 273.15;

Tprint = 210 + 273.15;
Tlow   = 65  + 273.15;
Thigh  = 205 + 273.15;

un_on = 24;
Rn    = 14.4;

Ln = 5e-2;
Cn = 1e-2;

Ct_n = 4.0;
Rt_n = 4.625;

% Initial nozzle condition
Tn  = Tprint;
Vc  = un_on;
iLn = un_on/Rn;

%% =========================================================
%  ELLIPSE
% ==========================================================

phi = @(t) omega*t;

xd = @(t) xc + a*cos(phi(t));
yd = @(t) yc + b*sin(phi(t));

dxd = @(t) -a*omega*sin(phi(t));
dyd = @(t)  b*omega*cos(phi(t));

ddxd = @(t) -a*omega^2*cos(phi(t));
ddyd = @(t) -b*omega^2*sin(phi(t));

dddxd = @(t)  a*omega^3*sin(phi(t));
dddyd = @(t) -b*omega^3*cos(phi(t));

% Required currents
ix = @(t) (rm/Km) .* ...
    (Meq*ddxd(t) + Bx*dxd(t));

iy = @(t) (rm/Km) .* ...
    (Meq*ddyd(t) + By*dyd(t));

% Current derivatives
dix = @(t) (rm/Km) .* ...
    (Meq*dddxd(t) + Bx*ddxd(t));

diy = @(t) (rm/Km) .* ...
    (Meq*dddyd(t) + By*ddyd(t));

% Required voltages
ux = @(t) Lm*dix(t) ...
    + Rm*ix(t) ...
    + (2*Km/rm)*dxd(t);

uy = @(t) Lm*diy(t) ...
    + Rm*iy(t) ...
    + (2*Km/rm)*dyd(t);

%% =========================================================
%  ELLIPSE LENGTH AND FILAMENT LENGTH
% ==========================================================

phi_grid = linspace(0,2*pi,10000);

ds_dphi = sqrt( ...
    a^2*sin(phi_grid).^2 + ...
    b^2*cos(phi_grid).^2);

s_grid = cumtrapz(phi_grid,ds_dphi);

Lellipse = s_grid(end);

% Every spool contains enough filament for 3/4 of an ellipse
Lspool = 0.75*Lellipse;

fprintf('Ellipse length: %.4f m\n',Lellipse);
fprintf('Filament length per spool: %.4f m\n\n',Lspool);

%% =========================================================
%  INITIAL CONDITIONS FOR HYBRID LOGIC
% ==========================================================

% We simulate two complete ellipses
phi_final = 4*pi;

phi_current = 0;

% Filament available in the first spool
Lremaining = Lspool;

% Logical variables
s   = 0;        % saved-position flag
ell = 1;        % filament loaded
p   = 1;        % deposition active

replacement = 0;

% Data used for the final figure
savedX = [];
savedY = [];

wipeX = {};
wipeY = {};

backX = {};
backY = {};

vpos = 0.10;

%% =========================================================
%  HYBRID SIMULATION
% ==========================================================

while phi_current < phi_final

    %% -----------------------------------------------------
    % q0 - PRINTING
    % ------------------------------------------------------

    % Position inside the current ellipse
    phi_local = mod(phi_current,2*pi);

    % Arc length already printed in the current ellipse
    s_current = interp1(phi_grid,s_grid,phi_local);

    % Distance from current point to the end of the ellipse
    LtoEnd = Lellipse - s_current;

    if Lremaining < LtoEnd

        % Filament finishes before completing the ellipse
        s_stop = s_current + Lremaining;

        phi_stop_local = interp1( ...
            s_grid,phi_grid,s_stop);

        currentLap = floor(phi_current/(2*pi));

        phi_stop = ...
            currentLap*2*pi + phi_stop_local;

        filamentFinished = true;

    else

        % Enough filament to complete the current ellipse
        currentLap = floor(phi_current/(2*pi));

        phi_stop = ...
            (currentLap + 1)*2*pi;

        filamentFinished = false;

    end

    % Do not go beyond the required number of ellipses
    phi_stop = min(phi_stop,phi_final);

    %% Simulate the printing motion

    t_start = phi_current/omega;
    t_stop  = phi_stop/omega;

    tPrint = linspace(t_start,t_stop,500)';

    U = [ux(tPrint),uy(tPrint)];

    X0 = [xd(t_start);
          dxd(t_start);
          ix(t_start);
          yd(t_start);
          dyd(t_start);
          iy(t_start)];

    Y=lsim(sys,U,tPrint-tPrint(1),X0);

    %% Update filament length

    if filamentFinished

        % All remaining filament has been used
        Lremaining = 0;

    else

        % Filament used to reach the end of the ellipse
        usedLength = LtoEnd;

        Lremaining = ...
            Lremaining - usedLength;
    end

    phi_current = phi_stop;

    %% Stop if the required two ellipses are complete

    if phi_current >= phi_final
        break
    end

    %% -----------------------------------------------------
    % FILAMENT FINISHED
    % -----------------------------------------------------

    if filamentFinished

        replacement = replacement + 1;

        %% q0 -> q1
        % Save current position

        xs = xd(phi_current/omega);
        ys = yd(phi_current/omega);

        s = 1;
        p = 0;

        savedX(replacement) = xs;
        savedY(replacement) = ys;

        fprintf('Replacement %d\n',replacement);
        fprintf('Saved position: x = %.3f m, y = %.3f m\n', ...
            xs,ys);

        %% -------------------------------------------------
        % q1 - STOP AND SAVE
        % -------------------------------------------------

        % The position has already been stored.
        % Since s = 1, the system can move to q2.

        %% -------------------------------------------------
        % q2 - MOVE TO WIPE TOWER
        % -------------------------------------------------

        [xWipe,yWipe] = moveHead( ...
            xs,ys,0,0,vpos, ...
            sys,rm,Km,Bx,By,Rm);

        wipeX{replacement} = xWipe;
        wipeY{replacement} = yWipe;

        %% -------------------------------------------------
        % q3 - WIPING
        % -------------------------------------------------

        p = 1;

        % Wiping finished: old filament has been removed
        ell = 0;

        p = 0;

        %% -------------------------------------------------
        % q4 - HEAT OFF
        % -------------------------------------------------

        [Tn,Vc,iLn] = changeNozzleTemperature( ...
            Tn,Vc,iLn, ...
            0,Tlow, ...
            Te,Rn,Ln,Cn,Ct_n,Rt_n, ...
            -1);

        %% -------------------------------------------------
        % q5 - LOAD NEW FILAMENT
        % -------------------------------------------------

        ell = 1;

        % A new spool is loaded
        Lremaining = Lspool;

        %% -------------------------------------------------
        % q6 - HEAT ON
        % -------------------------------------------------

        [Tn,Vc,iLn] = changeNozzleTemperature( ...
            Tn,Vc,iLn, ...
            un_on,Thigh, ...
            Te,Rn,Ln,Cn,Ct_n,Rt_n, ...
            1);

        %% -------------------------------------------------
        % q7 - MOVE BACK TO SAVED POSITION
        % -------------------------------------------------

        [xBack,yBack] = moveHead( ...
            0,0,xs,ys,vpos, ...
            sys,rm,Km,Bx,By,Rm);

        backX{replacement} = xBack;
        backY{replacement} = yBack;

        % The heater remains on during the return.
        % At the start of printing we assume that the nozzle
        % has reached its nominal printing condition again.
        Tn  = Tprint;
        Vc  = un_on;
        iLn = un_on/Rn;

        % Saved position reached
        s = 0;

        % Printing starts again
        p = 1;

        fprintf('Printing resumed.\n\n');

    end

end

%% =========================================================
%  PLOT OF THE TWO REPLACEMENTS
% ==========================================================

figure

phi_plot = linspace(0,2*pi,500);

for k = 1:min(2,replacement)

    subplot(1,2,k)
    hold on

    % Nominal ellipse
    plot( ...
        xc + a*cos(phi_plot), ...
        yc + b*sin(phi_plot), ...
        '--', ...
        'LineWidth',1.4)

    % Motion to wipe tower
    plot( ...
        wipeX{k}, ...
        wipeY{k}, ...
        'LineWidth',1.8)

    % Motion back to saved position
    plot( ...
        backX{k}, ...
        backY{k}, ...
        'LineWidth',1.8)

    % Position where printing was interrupted
    plot( ...
        savedX(k),savedY(k), ...
        'o', ...
        'MarkerSize',8, ...
        'LineWidth',1.6)

    % Wipe tower
    plot( ...
        0,0, ...
        's', ...
        'MarkerSize',8, ...
        'LineWidth',1.6)

    % Starting point of the ellipse
    plot( ...
        xc+a,yc, ...
        'ko', ...
        'MarkerFaceColor','k', ...
        'MarkerSize',5)

    axis equal
    grid on

    xlim([0 Lbed])
    ylim([0 Lbed])

    xlabel('x [m]')
    ylabel('y [m]')

    if k == 1
        title('First filament replacement')
    else
        title('Second filament replacement')
    end

end


%% =========================================================
%  LOCAL FUNCTIONS
% ==========================================================

function [xPath,yPath] = moveHead( ...
    x0,y0,xf,yf,vpos, ...
    sys,rm,Km,Bx,By,Rm)

    distance = sqrt( ...
        (xf-x0)^2 + ...
        (yf-y0)^2);

    if distance == 0
        xPath = x0;
        yPath = y0;
        return
    end

    Tmove = distance/vpos;

    vx = vpos*(xf-x0)/distance;
    vy = vpos*(yf-y0)/distance;

    ix = (rm/Km)*(Bx*vx);
    iy = (rm/Km)*(By*vy);

    ux = Rm*ix + (2*Km/rm)*vx;
    uy = Rm*iy + (2*Km/rm)*vy;

    t = linspace(0,Tmove,300)';

    U = [ ...
        ux*ones(size(t)), ...
        uy*ones(size(t))];

    X0 = [ ...
        x0;
        vx;
        ix;
        y0;
        vy;
        iy];

    Y = lsim(sys,U,t,X0);

    xPath = Y(:,1);
    yPath = Y(:,2);

end


function [T,Vc,iL] = changeNozzleTemperature( ...
    T0,Vc0,iL0, ...
    un,Ttarget, ...
    Te,Rn,Ln,Cn,Ct,Rt,direction)

    z0 = [T0;Vc0;iL0];

    options = odeset( ...
        'Events', ...
        @(t,z) temperatureEvent( ...
        t,z,Ttarget,direction));

    [~,z] = ode45( ...
        @(t,z) nozzleDynamics( ...
        t,z,un,Te,Rn,Ln,Cn,Ct,Rt), ...
        [0 500], ...
        z0, ...
        options);

    T  = z(end,1);
    Vc = z(end,2);
    iL = z(end,3);

end

function dz = nozzleDynamics( ...
    ~,z,un,Te,Rn,Ln,Cn,Ct,Rt)

    T  = z(1);
    Vc = z(2);
    iL = z(3);

    dT = ...
        Vc^2/(Rn*Ct) ...
        - (T-Te)/(Ct*Rt);

    dVc = ...
        iL/Cn ...
        - Vc/(Rn*Cn);

    diL = ...
        -Vc/Ln ...
        + un/Ln;

    dz = [dT;dVc;diL];

end

function [value,isterminal,direction] = ...
    temperatureEvent(~,z,Ttarget,eventDirection)

    value = z(1) - Ttarget;

    isterminal = 1;
    direction = eventDirection;

end