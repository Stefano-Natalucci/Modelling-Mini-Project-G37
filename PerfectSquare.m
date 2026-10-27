clear; clc; close all;

%% Parameters

rm = 0.006;              % Motor pulley radius [m]
Ls = 0.20;               % Square side length [m]
v_print = 0.01;          % Printing speed [m/s]
dt = 0.01;               % Time step [s]

theta_s = Ls/rm;         % Required angular displacement for one side
omega = v_print/rm;      % Motor angular speed [rad/s]

% Initial position of the square
x0 = 0.1;
y0 = 0.1;

% Initial motor angles from CoreXY kinematics
theta_l = (x0 + y0)/rm;
theta_r = (x0 - y0)/rm;

% Reference angles at the beginning of the first side
theta_l_s = theta_l;
theta_r_s = theta_r;

%% Hybrid modes
% q = 1 : qR  -> H+
% q = 2 : qU  -> V+
% q = 3 : qL  -> H-
% q = 4 : qD  -> V-
% q = 5 : qF  -> Square completed

q = 1;

%% Save simulation data

theta_l_hist = theta_l;
theta_r_hist = theta_r;

%% Simulation

while q ~= 5

    % Motor angular velocities according to the active mode
    switch q

        case 1      % qR - first side, H+
            omega_l =  omega;
            omega_r =  omega;

        case 2      % qU - second side, V+
            omega_l =  omega;
            omega_r = -omega;

        case 3      % qL - third side, H-
            omega_l = -omega;
            omega_r = -omega;

        case 4      % qD - fourth side, V-
            omega_l = -omega;
            omega_r =  omega;
    end

    % Continuous evolution
    theta_l = theta_l + omega_l*dt;
    theta_r = theta_r + omega_r*dt;

    % Angular displacement measured from the beginning of the side
    Delta_theta_l = theta_l - theta_l_s;
    Delta_theta_r = theta_r - theta_r_s;

    %% Guards and resets

    switch q

        case 1
            % qR -> qU
            if Delta_theta_l >= theta_s && Delta_theta_r >= theta_s
                q = 2;

                theta_l_s = theta_l;
                theta_r_s = theta_r;
            end

        case 2
            % qU -> qL
            if Delta_theta_l >= theta_s && Delta_theta_r <= -theta_s
                q = 3;

                theta_l_s = theta_l;
                theta_r_s = theta_r;
            end

        case 3
            % qL -> qD
            if Delta_theta_l <= -theta_s && Delta_theta_r <= -theta_s
                q = 4;

                theta_l_s = theta_l;
                theta_r_s = theta_r;
            end

        case 4
            % qD -> qF
            if Delta_theta_l <= -theta_s && Delta_theta_r >= theta_s
                q = 5;
            end
    end

    % Save angular positions
    theta_l_hist(end+1) = theta_l;
    theta_r_hist(end+1) = theta_r;

end

%% Reconstruct Cartesian position only for visualisation

x = rm/2 * (theta_l_hist + theta_r_hist);
y = rm/2 * (theta_l_hist - theta_r_hist);

%% Plot

figure;

plot(x, y, 'LineWidth', 2);
hold on;

plot(x(1), y(1), 'o', 'MarkerSize', 7);
plot(x(end), y(end), 's', 'MarkerSize', 7);

xlabel('x [m]');
ylabel('y [m]');
title('Square Drawing Using Motor Angular Measurements');

legend('Printed square', 'Starting point', 'Final point', ...
       'Location', 'best');

axis equal;
grid on;

xlim([0 0.40]);
ylim([0 0.40]);