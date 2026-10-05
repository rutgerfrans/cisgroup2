%% Physical parameters used by Simscape
%%%%% You can vary the parameters %%%%%%
param.m  = 1.0;     % Pendulum mass [kg]
param.M  = 1.0;     % Cart mass [kg]
param.L  = 0.5;     % Full pole length [m]
param.bc = 0.1;     % Cart damping [N*s/m]
param.bp = 0.05;    % Pivot damping [N*m*s/rad]
g=9.81;

%% Initial conditions: hanging down, stationary
%%%%% Vary the initial state %%%%%%
param.p0     = 0;
param.v0     = 0;
param.theta0 = 10; % Simscape joint target is in degrees (0 upright, 180 down)
param.omega0 = 0;

x0 = [param.p0; param.v0;
      deg2rad(param.theta0);param.omega0];

%% Parameters for nonlinear dynamics 
model.M   = param.M;
model.m   = param.m;
model.ell = param.L/2;     % Pivot-to-point-mass distance [m]
model.bc  = param.bc;
model.bp  = param.bp;
model.g   = g;

%%%%% You can vary these parameters %%%%%%
Ts   = 0.01;               % Sample time for discrete-time dynamics [s]
Fmax = 20;                 % Force saturation [N]
xmax = 0.5;                % Cart-travel limit [m] (can be arbitrarily large)

%% Continuous-time linearization at upright
% Upright = 0; positive angles are counterclockwise from 0
%%% You do not need to change this part %%%

M   = model.M;
m   = model.m;
ell = model.ell;

J = m*ell^2;
Delta = (M+m)*J - (m*ell)^2;

A = [0, 1, 0, 0;
     0, -J*model.bc/Delta, ...
         m^2*g*ell^2/Delta, -m*ell*model.bp/Delta;
     0, 0, 0, 1;
     0, -m*ell*model.bc/Delta, ...
         (M+m)*m*g*ell/Delta, -(M+m)*model.bp/Delta];

B = [0;
     J/Delta;
     0;
     m*ell/Delta];

%% Discrete-time LQR (Task 2.1)

C = eye(4,4);
D = zeros(4,1);

sys_cont = ss(A,B,C,D);
sys_disc = c2d(sys_cont, Ts, 'zoh');

%%% Obtain a discrete-time model of the continuous time system
Ad = sys_disc.A;
Bd = sys_disc.B;

% ---- 2.1.1 confirmation ---- %
% these should equal
% sort(eig(Ad))
% sort(exp(eig(A)*Ts))


%%% Compute the discrete LQR gain, select and justify weight matrices

% LQR weights: Bryson's rule weight=1/(maximum acceptable deviation^2)

v_max = 1; % 1 meter per second 
theta_max_deg = 5; % small degree since model is linearized upright
theta_max = deg2rad(theta_max_deg); % convert in radians
omega_max = 1; % recovery speed from small tilts

q_p = 1 / xmax^2; % cart position weight in m
q_v = 1 / v_max^2; % cart velocity weight in m/s
q_theta = 1 / theta_max^2; % pole angle weight
q_omega = 1 / omega_max^2; % angular velocity weight


Q_lqr = diag([q_p q_v q_theta q_omega]);
R_lqr = 1/(Fmax^2);

K_lqr = dlqr(Ad, Bd, Q_lqr, R_lqr);

% ---- 2.1.2 confirmation ---- %
% size(K_lqr) % should be 1 x 4
% abs(eig(Ad - Bd*K_lqr)) % should be stable -> all <1

% temporary values for testing upright only
% should later be replaced by iLQR output values
Knom = 0;
Unom = 0;
Xnom = 0;


%% iLQR setup
% T = ...           % Swing-up horizon [s]
% N = round(T/Ts);
% T = N*Ts;

% Initial guess
% U0 = ...

% Swing-up weight matrices
% Q  = ...
% R  = ...
% Qf = ...

%% Optimize swing-up
%%% Implement ilqr function %%%
[Xnom,Unom,Knom,costHistory] = ilqr(x0,U0,model,Ts,Q,R,Qf,Fmax);

%% Check nominal trajectory

%%% Check if the nominal trajectory Xnum, Unom returned by iLQR satisfies
%%% that the terminal state is in the capture zone, the cart displacement
%%% is within limits, the force is bounded (it should be)


%% Plot results

%%%% You can plot the nominal state and control trajectories. These can be
%%%% informative for changing iLQR parameters for finding a feasible
%%%% trajectory.

%% Simulation 
out=sim("cartpole");

%%% Plot additional results. For example, comparing nominal and resulting
%%% swing-up trajectory.

