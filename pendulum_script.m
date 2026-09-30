%% Physical parameters 
%%%%% These are the parameters used by the Simulink model and they
%%%%% correspond to the properties of the real physical system

p.m = 1.0;        % kg
p.l = 0.5;        % m, distance pivot -> point mass
p.b = 0.05;       % N*m/(rad/s)
g = 9.81;       % m/s^2

p.theta0 = 35;  %degrees
p.omega0 = 0;   %rad/s

%% Model parameters 
%%%%%%% These are the parameters you can use for your model. Currently,
%%%%%%% they are set equal to the values of the physical system but you can
%%%%%%% change and observe how your controller and observer behave under model
%%%%%%% uncertainty

m.m = p.m;        % kg
m.l = p.l;        % m, distance pivot -> point mass
m.b = p.b;       % N*m/(rad/s)
m.J = m.m*m.l^2;


%% Sampling
%%%%%% Controller works with measurements obtained at regular intervals.
%%%%%% The control inputs are in discrete time as well. They are then
%%%%%% applied to the system in zero-order-hold manner. 
%%%%%% Sampling period: every Ts seconds
Ts = 0.02; 

%% 1.1
u = 0;
figure; hold on; grid on
for th0 = [0 0.01 -0.01]
    [t,X] = ode45(@(t,x) pen_sys(t,x,p,g,u), [0 50], [th0; 0]);
    plot(t, rad2deg(X(:,1)))
end
xlabel('t [s]'); ylabel('\theta [deg]')

%% Actuator limits 
umin = -2;
umax = 2;

%% Linear model (1.2)
%%%%%% Here you should set the matrices Ac and Bc which correspond to
%%%%%% the linearized system 

Ac = [0 1; (m.m * g * m.l)/m.J -m.b/m.J];
Bc = [0; 1/m.J];

eigenvalues_of_Ac = eig(Ac);
disp("Eigenvalues of Ac:");
disp(eigenvalues_of_Ac);
% eigenvalues are 4.3306 and -4.5306
% one has a positive real part, so the system is unstable
% this is similar to the unstability of the pendulum when the initial theta
% is -0.01 or 0.01

%% Discrete-time model (1.3)
%%%%% Disretize the linear apprroximation; you can use MATLAB functions
%%%%% from control systtem toolbox such as c2d or compute matrix
%%%%% exponentials

sys_c = ss(Ac, Bc, eye(2), zeros(2,1));

sys_d = c2d(sys_c, Ts, 'zoh');

Ad = sys_d.A;
Bd = sys_d.B;

%% LQR
%%%%% Design LQR 

Q = [1 0; 0 1];
R = 1;

[K, S, P] = dlqr(Ad, Bd, Q, R);

%% Observer
%%%% Design an observer using pole placement
C = [1 0];
%xhat0 = ... 

%L = ...

%% Simulate the actuated pendulum
%%% The following command runs the simulation.
%%% You should see an animation of the pendulum and out is a struct that
%%% contains the outputs from the simulation (for example, measured and
%%% estimated values, control input, or any other thing you add)

out=sim("actuated_pendulum");

theta_sim = out.theta_sim;
omega_sim = out.omega_sim;
t = theta_sim.Time;

figure;

ax1 = subplot(2, 1, 1);
plot(t, theta_sim.Data)
title("Simulated Theta");
xlabel('t [s]'); ylabel('\theta [deg]')
grid on

ax2 = subplot(2, 1, 2);
plot(t, omega_sim.Data)
title("Simulated Omega");
xlabel('t [s]'); ylabel('\omega [deg/s]')
grid on

sgtitle("Simulation outputs with noise at theta0=" + string(p.theta0) + " omega0=" + ...
    string(p.omega0) + " umin=" + string(umin) + " umax=" + string(umax) + ...
    newline + "Q=" + string(mat2str(Q,3)) + " R=" + string(R) + " length=2m")

%% First order equation
function xdot = pen_sys(t, x, p, g, u)
theta = x(1);
omega = x(2);

J = p.m*p.l^2;
theta_ddot = (p.m * g * p.l * sin(theta) - p.b * omega + u) / J;

xdot = [omega;theta_ddot];
end


