clear
close all
clc

a = 0.10;      % semi-major axis [m]
b = 0.05;      % semi-minor axis [m]
T = 10;        % period [s]
omega = 2*pi/T;

xd   = @(t) a*cos(omega*t);
yd   = @(t) b*sin(omega*t);

dxd  = @(t) -a*omega*sin(omega*t);
dyd  = @(t)  b*omega*cos(omega*t);

ddxd = @(t) -a*omega^2*cos(omega*t);
ddyd = @(t) -b*omega^2*sin(omega*t);

% Parameters
m     = 0.01;
Jm    = 0.00005;
rm    = 0.02;
Km    = 1;

mu_r  = 0.001;
mu_lx = 0.005;
mu_ly = 0.01;

Lm = 1e-3;
Rm = 100;

% Equivalent parameters
Meq = m + 2*Jm/rm^2;

Bx = 2*mu_lx + 8*mu_r;
By = 2*mu_ly + 4*mu_r;

% Required transformed currents
ix = @(t) (rm/Km) * ...
    (Meq*ddxd(t) + Bx*dxd(t));

iy = @(t) (rm/Km) * ...
    (Meq*ddyd(t) + By*dyd(t));

% Motor currents
iL = @(t) 0.5*(ix(t) + iy(t));
iR = @(t) 0.5*(ix(t) - iy(t));

% Third derivatives of desired trajectory
dddxd = @(t)  a*omega^3*sin(omega*t);
dddyd = @(t) -b*omega^3*cos(omega*t);

% Current derivatives
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

% Motor voltages
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

plot(x,y)
axis equal
grid on
xlabel('x [m]')
ylabel('y [m]')