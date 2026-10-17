clear
close all
clc

a = 0.08;          % scale [m]
T = 10;            % period [s]
omega = 2*pi/T;

xc = 0;            % set e.g. 0.145 to translate the shape
yc = 0;

phi = @(t) omega*t;

% Desired trajectory
xd = @(t) xc + a*sin(2*phi(t)).*cos(phi(t));
yd = @(t) yc + a*sin(2*phi(t)).*sin(phi(t));

% First derivatives
dxd = @(t) a*omega*( ...
      2*cos(2*phi(t)).*cos(phi(t)) ...
    - sin(2*phi(t)).*sin(phi(t)) );

dyd = @(t) a*omega*( ...
      2*cos(2*phi(t)).*sin(phi(t)) ...
    + sin(2*phi(t)).*cos(phi(t)) );

% Second derivatives
ddxd = @(t) a*omega^2*( ...
     -5*sin(2*phi(t)).*cos(phi(t)) ...
     -4*cos(2*phi(t)).*sin(phi(t)) );

ddyd = @(t) a*omega^2*( ...
     -5*sin(2*phi(t)).*sin(phi(t)) ...
     +4*cos(2*phi(t)).*cos(phi(t)) );

% Third derivatives
dddxd = @(t) a*omega^3*( ...
     -14*cos(2*phi(t)).*cos(phi(t)) ...
     +13*sin(2*phi(t)).*sin(phi(t)) );

dddyd = @(t) a*omega^3*( ...
     -14*cos(2*phi(t)).*sin(phi(t)) ...
     -13*sin(2*phi(t)).*cos(phi(t)) );

% Physical parameters
m     = 0.01;
Jm    = 0.00005;
rm    = 0.02;
Km    = 1;

mu_r  = 0.001;
mu_lx = 0.005;
mu_ly = 0.01;

Lm = 1e-3;      % motor inductance [H]
Rm = 100;       % motor resistance [Ohm]

% Equivalent parameters
Meq = m + 2*Jm/rm^2;

Bx = 2*mu_lx + 8*mu_r;
By = 2*mu_ly + 4*mu_r;

% Required transformed currents
ix = @(t) (rm/Km) * ...
    (Meq*ddxd(t) + Bx*dxd(t));

iy = @(t) (rm/Km) * ...
    (Meq*ddyd(t) + By*dyd(t));

% Physical motor currents
iL = @(t) 0.5*(ix(t) + iy(t));
iR = @(t) 0.5*(ix(t) - iy(t));

% Derivatives of transformed currents
dix = @(t) (rm/Km) * ...
    (Meq*dddxd(t) + Bx*ddxd(t));

diy = @(t) (rm/Km) * ...
    (Meq*dddyd(t) + By*ddyd(t));

% Required transformed voltages
ux = @(t) Lm*dix(t) ...
    + Rm*ix(t) ...
    + (2*Km/rm)*dxd(t);

uy = @(t) Lm*diy(t) ...
    + Rm*iy(t) ...
    + (2*Km/rm)*dyd(t);

% Physical motor voltages
uL = @(t) 0.5*(ux(t) + uy(t));
uR = @(t) 0.5*(ux(t) - uy(t));

% Decoupled state-space model
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

C = [1 0 0 0 0 0;
     0 0 0 1 0 0];

D = zeros(2,2);

sys = ss(A,B,C,D);

% Simulation
t = linspace(0,T,1000)';

Ux = ux(t);
Uy = uy(t);

X0 = [xd(0);
      dxd(0);
      ix(0);
      yd(0);
      dyd(0);
      iy(0)];

U = [Ux Uy];

[Y,t,X] = lsim(sys,U,t,X0);

x = Y(:,1);
y = Y(:,2);

% Plot
plot(x,y)
axis equal
grid on
xlabel('x [m]')
ylabel('y [m]')