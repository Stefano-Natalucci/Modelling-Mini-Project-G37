%% Defining the Parameters for the model

clear;
clc;
close all;

w_fil = 0.0004;        % [m] filament side length
rho = 1300;            % [kg/m^3]
C = 1800;              % [J/(kg*K)]
k_filament = 0.13;     % [W/(m*K)]

Tn = 210 + 273.15;     % [K] nozzle temperature
Tb = 60 + 273.15;      % [K] heated bed temperature
Te = 25 + 273.15;      % [K] ambient temperature

h_air = 20;            % [W/(m^2*K)]
h_plate = 100;         % [W/(m^2*K)]

dx = 0.0001;           % [m]
dt = 0.01;             % [s] (dt has been chosen keeping in mind the condition for stability in PDEs)

L_tot = 0.05;          % [m]
v_print = 0.01;        % [m/s]

tolerance = 5;         % [K]

N_tot = round(L_tot/dx);  %Total number of printed nodes
t_print = L_tot/v_print;   %Total printing time

alpha = k_filament/(rho*C);  

% Creating the initiLal T vector onctaining only room temperature values.
T = Te * ones(1,N_tot);

%% PDE While Printing

t = 0;

while t < t_print

    % Number of deposited nodes at the current time
    N_done = min(floor(v_print*t/dx) + 1, N_tot);

    % Temperature field used for the current update
    T_old = T;

    % Newly deposited filament is at nozzle temperature
    T_old(N_done) = Tn;

    T_temp = T_old;

    % If at least two nodes have been deposited,
    % update the left boundary using zero heat-flux condition
    if N_done >= 2

        % Neumann condition dT/dx = 0 at x = 0
        % implemented with a ghost-node approximation
        T_left_ghost = T_old(2);

        T_temp(1) = T_old(1) + dt * ( ...
            alpha * (T_old(2) - 2*T_old(1) + T_left_ghost)/(dx^2) ...
            + ((3*w_fil)*h_air*(Te - T_old(1)) ...
            + w_fil*h_plate*(Tb - T_old(1))) ...
            /(rho*C*w_fil^2) );

    end

    % Interior deposited nodes
    for i = 2:(N_done-1)

        T_temp(i) = T_old(i) + dt * ( ...
            alpha * (T_old(i+1) - 2*T_old(i) + T_old(i-1))/(dx^2) ...
            + ((3*w_fil)*h_air*(Te - T_old(i)) ...
            + w_fil*h_plate*(Tb - T_old(i))) ...
            /(rho*C*w_fil^2) );

    end

    % The active printing front remains at nozzle temperature
    T_temp(N_done) = Tn;

    T = T_temp;
    t = t + dt;

end

T_hot = T;

fprintf('Printing time: %.2f s\n', t_print);

%% PDE in Closed Loop (Printing is Finished)

T_max = max(T);
T_min = min(T);

t_cooling = 0;

while (T_max - T_min > tolerance)

    T_temp = T;

    for i = 1:N_tot

        % Periodic left neighbour
        i_left = i - 1;
        if i_left < 1
            i_left = N_tot;
        end

        % Periodic right neighbour
        i_right = i + 1;
        if i_right > N_tot
            i_right = 1;
        end

        T_temp(i) = T(i) + dt * (alpha * (T(i_right) - 2*T(i) + T(i_left))/(dx^2) ...
            + ((3*w_fil)*h_air*(Te - T(i)) ...
            + w_fil*h_plate*(Tb - T(i)))/(rho*C*w_fil^2) );

    end

    T = T_temp;

    T_max = max(T);
    T_min = min(T);

    t = t + dt;
    t_cooling = t_cooling + dt;

end

T_cold = T;

fprintf('Cooling time required to obtain a uniform temperature with tolerance %.2f K: %.3f s\n', ...
         tolerance, t_cooling);

fprintf('Total simulation time: %.3f s\n', t);

%% Temperature Comparison Between Hot and Cold Filament

figure('Name','Temperature Profile Comparison','Color','w');

x_axis = (0:N_tot-1) * dx * 1000;

plot(x_axis, T_hot - 273.15, 'LineWidth', 2, ...
    'DisplayName', sprintf('End of printing (t = %.1f s)', t_print));

hold on;

plot(x_axis, T_cold - 273.15, 'LineWidth', 2, ...
    'DisplayName', sprintf('Nearly uniform (t = %.1f s)', t));

yline(Tb - 273.15, '--', 'Bed temperature', ...
    'LabelHorizontalAlignment','left', ...
    'HandleVisibility','off');

yline(Te - 273.15, ':', 'Ambient temperature', ...
    'LabelHorizontalAlignment','left', ...
    'HandleVisibility','off');

yline(Tn - 273.15, ':', 'Nozzle temperature', ...
    'LabelHorizontalAlignment','left', ...
    'HandleVisibility','off');

grid on;

xlabel('Position along the filament [mm]');
ylabel('Temperature [^\circC]');
title('Filament Temperature at the End of Printing and After Cooling');

legend('Location','northeast');