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
param.theta0 = 150; % Simscape joint target is in degrees (0 upright, 180 down)
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
Fmax = 40;                 % Force saturation [N]
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

%% Empirical Capture Zone - angle vs. angular velocity (Task 2.2.1)

thetas = deg2rad(-90:1:90); % list of all theta values spaced by 1 degree
omegas = -10:0.25:10; % list of all omega values spaced by 0.25
Tsim = 10; % run each sim 10 seconds

captured_map = false(length(omegas), length(thetas)); % initialize all theta-omega pairings as false
for i = 1 : length(omegas) 
    for j = 1 : length(thetas) % nested loop through both values
        x0_test = [0; 0; thetas(j); omegas(i)]; % initialize test conditions x0
        captured_map(i, j) = test_capture(x0_test, K_lqr, model, Ts, Fmax, xmax, Tsim); % run capture zone experiment
    end
end

% visualize capture zone
figure; imagesc(rad2deg(thetas), omegas, captured_map)
set(gca,'YDir','normal'); colormap([1 0.6 0.6; 0.6 0.9 0.6])
xlabel('\theta_0 [deg]'); ylabel('\omega_0 [rad/s]')
title('Capture Zone: green = stabilized, red = unstable, p_0 = v_0 = 0')

%% Empirical Capture Zone - angle vs. cart position (Task 2.2.2)
p_list  = -xmax:0.05:xmax; % list of all cart velocities spaced by 0.05

captured_p = false(length(p_list), length(thetas));
for i = 1:length(p_list)
    for j = 1:length(thetas)
        x0_test = [p_list(i); 0; thetas(j); 0];
        captured_p(i,j) = test_capture(x0_test, K_lqr, model, Ts, Fmax, xmax, Tsim);
    end
end

figure; imagesc(rad2deg(thetas), p_list, double(captured_p))
set(gca,'YDir','normal'); colormap([1 0.6 0.6; 0.6 0.9 0.6])
xlabel('\theta_0 [deg]'); ylabel('p_0 [m]')
title('Capture zone (green = stabilized), v_0 = \omega_0 = 0')

%% Empirical Capture Zone - angle vs. cart velocity (Task 2.2.3)
v_list = -4:0.25:4; % list of all cart velocities spaced by 0.25

captured_v = false(length(v_list), length(thetas));
for i = 1:length(v_list)
    for j = 1:length(thetas)
        x0_test = [0; v_list(i); thetas(j); 0];
        captured_v(i,j) = test_capture(x0_test, K_lqr, model, Ts, Fmax, xmax, Tsim);
    end
end

figure; imagesc(rad2deg(thetas), v_list, double(captured_v))
set(gca,'YDir','normal'); colormap([1 0.6 0.6; 0.6 0.9 0.6])
xlabel('\theta_0 [deg]'); ylabel('v_0 [m/s]')
title('Capture zone (green = stabilized), p_0 = \omega_0 = 0')

%% Capture Criterion (Taks 2.2.4)
% manually read off limits from earlier capture zone plots to determine
% rough limits --> later toyed with to arrive at these capture zone limits
cap.p = 0.11;
cap.v = 0.3;
cap.theta = deg2rad(10);
cap.omega = 1.11;
capture_limits = [cap.p; cap.v; cap.theta; cap.omega]; % limit box

% test every combination at values -limit, 0, +limit (middle + extreme points)
[a,b,c,d] = ndgrid([-1 0 1]);
S = [a(:) b(:) c(:) d(:)]';

corner_check = false(1, size(S,2)); % initially false matrix for all combinations
for k = 1 : size(S, 2) % loop through all possible values
    x0_test = S(:,k) .* capture_limits; % apply box limits for this iteration
    corner_check(k) = test_capture(x0_test, K_lqr, model, Ts, Fmax, xmax, Tsim); % test capture, store in corner_check
end
fprintf('Captured: %d of %d\n', sum(corner_check), numel(corner_check)); % % of combinations succeeded

disp(S(:, ~corner_check)') % print failed combinations, used for manual tweaking of limits

%% iLQR setup
T = 8;           % Swing-up horizon [s]
N = round(T/Ts);
T = N*Ts;

% Initial guess
U0 = 50 * ones(1,N);

% Swing-up weight matrices
Q  = [0.4, 0.2, 10, 1] .* eye(4);
R  = 0.001;
Qf = 100000 * [4, 0.2, 10, 4] .* eye(4);

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

