clear
close all
clc

%% =========================================================
%  PRINTING AREA AND CLOVER GEOMETRY
% ==========================================================

Lbed = 0.40;          % print-bed side length [m]

xc = Lbed/2;
yc = Lbed/2;

a = 0.15;             % clover scale [m]

omega = 2*pi/10;      % angular speed [rad/s]

% Starting point of the clover (phi = 0)
xs = xc;
ys = yc;

%% =========================================================
%  PHYSICAL PARAMETERS
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

% Equivalent parameters
Meq = m + 2*Jm/rm^2;

Bx = 2*mu_lx + 8*mu_r;
By = 2*mu_ly + 4*mu_r;

%% =========================================================
%  STATE-SPACE MODEL
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

C = [1 0 0 0 0 0;
     0 0 0 1 0 0];

D = zeros(2,2);

sys = ss(A,B,C,D);

%% =========================================================
%  PHASE 1 - POSITIONING TO START POINT
% ==========================================================

vpos = 0.10;                  % constant positioning speed [m/s]

distance_out = sqrt(xs^2 + ys^2);
Tmove = distance_out/vpos;

vx_out = vpos*xs/distance_out;
vy_out = vpos*ys/distance_out;

% Constant transformed currents
ix_move = (rm/Km)*(Bx*vx_out);
iy_move = (rm/Km)*(By*vy_out);

% Constant transformed voltages
ux_move = Rm*ix_move + (2*Km/rm)*vx_out;
uy_move = Rm*iy_move + (2*Km/rm)*vy_out;

% Physical motor quantities
iL_move = 0.5*(ix_move + iy_move);
iR_move = 0.5*(ix_move - iy_move);

uL_move = 0.5*(ux_move + uy_move);
uR_move = 0.5*(ux_move - uy_move);

t1 = linspace(0,Tmove,400)';

U1 = [ux_move*ones(size(t1)), ...
      uy_move*ones(size(t1))];

X01 = [0;
       vx_out;
       ix_move;
       0;
       vy_out;
       iy_move];

[Y1,t1,X1] = lsim(sys,U1,t1,X01);

x1 = Y1(:,1);
y1 = Y1(:,2);

%% =========================================================
%  PHASE 2 - PRE-PRINT WAIT
% ==========================================================

Twait_pre = 1;                % waiting time before printing [s]

t2 = linspace(0,Twait_pre,100)';

U2 = zeros(length(t2),2);

X02 = [xs;
       0;
       0;
       ys;
       0;
       0];

[Y2,t2,X2] = lsim(sys,U2,t2,X02);

x2 = Y2(:,1);
y2 = Y2(:,2);

%% =========================================================
%  PHASE 3 - PRINT CLOVER
% ==========================================================

phi_end = 2*pi;
Tclover = phi_end/omega;

phi = @(t) omega*t;

% Desired clover trajectory
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

% Required transformed currents
ix = @(t) (rm/Km) * ...
    (Meq*ddxd(t) + Bx*dxd(t));

iy = @(t) (rm/Km) * ...
    (Meq*ddyd(t) + By*dyd(t));

% Physical motor currents
iL = @(t) 0.5*(ix(t) + iy(t));
iR = @(t) 0.5*(ix(t) - iy(t));

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

% Physical motor voltages
uL = @(t) 0.5*(ux(t) + uy(t));
uR = @(t) 0.5*(ux(t) - uy(t));

t3 = linspace(0,Tclover,1000)';

U3 = [ux(t3), uy(t3)];

X03 = [xd(0);
       dxd(0);
       ix(0);
       yd(0);
       dyd(0);
       iy(0)];

[Y3,t3,X3] = lsim(sys,U3,t3,X03);

x3 = Y3(:,1);
y3 = Y3(:,2);

%% =========================================================
%  PHASE 4 - POST-PRINT WAIT
% ==========================================================

Twait_post = 1;               % waiting time after printing [s]

xf = x3(end);
yf = y3(end);

t4 = linspace(0,Twait_post,100)';

U4 = zeros(length(t4),2);

X04 = [xf;
       0;
       0;
       yf;
       0;
       0];

[Y4,t4,X4] = lsim(sys,U4,t4,X04);

x4 = Y4(:,1);
y4 = Y4(:,2);

%% =========================================================
%  PHASE 5 - RETURN TO ORIGIN
% ==========================================================

distance_home = sqrt(xf^2 + yf^2);
Thome = distance_home/vpos;

vx_home = -vpos*xf/distance_home;
vy_home = -vpos*yf/distance_home;

% Constant transformed currents
ix_home = (rm/Km)*(Bx*vx_home);
iy_home = (rm/Km)*(By*vy_home);

% Constant transformed voltages
ux_home = Rm*ix_home + (2*Km/rm)*vx_home;
uy_home = Rm*iy_home + (2*Km/rm)*vy_home;

% Physical motor quantities
iL_home = 0.5*(ix_home + iy_home);
iR_home = 0.5*(ix_home - iy_home);

uL_home = 0.5*(ux_home + uy_home);
uR_home = 0.5*(ux_home - uy_home);

t5 = linspace(0,Thome,400)';

U5 = [ux_home*ones(size(t5)), ...
      uy_home*ones(size(t5))];

X05 = [xf;
       vx_home;
       ix_home;
       yf;
       vy_home;
       iy_home];

[Y5,t5,X5] = lsim(sys,U5,t5,X05);

x5 = Y5(:,1);
y5 = Y5(:,2);

%% =========================================================
%  PRINTING BOOLEAN
% ==========================================================

printing1 = zeros(size(t1));
printing2 = zeros(size(t2));
printing3 = ones(size(t3));
printing4 = zeros(size(t4));
printing5 = zeros(size(t5));

printing = [printing1;
            printing2;
            printing3;
            printing4;
            printing5];

%% =========================================================
%  GLOBAL TIME VECTOR
% ==========================================================

t1_global = t1;

t2_global = Tmove + t2;

t3_global = Tmove ...
          + Twait_pre ...
          + t3;

t4_global = Tmove ...
          + Twait_pre ...
          + Tclover ...
          + t4;

t5_global = Tmove ...
          + Twait_pre ...
          + Tclover ...
          + Twait_post ...
          + t5;

t_global = [t1_global;
            t2_global;
            t3_global;
            t4_global;
            t5_global];

%% =========================================================
%  COMPLETE TRAJECTORY
% ==========================================================

x_total = [x1;
           x2;
           x3;
           x4;
           x5];

y_total = [y1;
           y2;
           y3;
           y4;
           y5];

%% =========================================================
%  PLOT
% ==========================================================

figure
hold on

% Print-bed boundary
plot([0 Lbed Lbed 0 0], ...
     [0 0 Lbed Lbed 0], ...
     'k-', 'LineWidth', 1)

% Positioning motion
plot(x1,y1,'--','LineWidth',1.5)

% Start of printing
plot(xs,ys,'o','MarkerSize',7,'LineWidth',1.5)

% Printed clover
plot(x3,y3,'-','LineWidth',2)

% Return to origin
plot(x5,y5,'--','LineWidth',1.5)

% Origin
plot(0,0,'s','MarkerSize',7,'LineWidth',1.5)

axis equal
grid on

xlim([0 Lbed])
ylim([0 Lbed])

xlabel('x [m]')
ylabel('y [m]')

legend('Print bed', ...
       'Positioning motion', ...
       'Start of printing', ...
       'Printed clover', ...
       'Return to origin', ...
       'Origin', ...
       'Location','best')

title('Simulated printing of a clover trajectory')