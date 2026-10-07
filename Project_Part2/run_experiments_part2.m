%% Experiments part 2 (task 4.1 and 5.1)
% every iLQR solution is also simulated on the Simscape model (cartpole.slx)

clear; close all;
out = '../experiment_results/part2/';
mkdir([out 'figures']);

% parameters and K_lqr from the first part of setup_cartpole.m
% (the capture zone sweeps after that take too long to rerun every time)
src = fileread('setup_cartpole.m');
idx = strfind(src, '%% Empirical Capture Zone');
eval(src(1:idx(1)-1));
x0 = [0; 0; pi; 0]; % hanging down
cap = [0.11; 0.3; deg2rad(10); 1.11]; % capture zone from task 2.2
[Ad_down, Bd_down] = linearize_trajectory([0; 0; pi; 0], 0, model, Ts);
Kdown = dlqr(Ad_down, Bd_down, diag([10 1 10 1]), 0.01); % for task 5.1
dF_amp = 0; dF_t0 = 6; dF_dur = 0.1; % no push on the cart
wrap = @(a) atan2(sin(a), cos(a));

% baseline settings
base.T = 3; base.Fmax = 40;
base.Q = diag([5 0.1 1 0.1]); base.R = 0.01; base.Qf = diag([1e4 1e3 1e4 1e3]);
base.U0 = @(t) 0*t; base.x0 = x0;
Qf = base.Qf;

%% list of experiments: {group, name, {field, value, ...}}
exps = {};
% 4.1 experiment 1: varying the problem, maximum force
for F = [10 20 25 30 35 40 60 80]
    exps(end+1,:) = {'E1a_max_force', sprintf('Fmax = %d N', F), {'Fmax', F}};
end
% 4.1 experiment 1: varying the problem, allowed cart displacement (via q_p)
for qp = [0 1 5 20 100]
    exps(end+1,:) = {'E1b_cart_displacement', sprintf('q_p = %d', qp), {'Q', diag([qp 0.1 1 0.1])}};
end
% 4.1 experiment 2: choosing Q, Qf and R
exps = [exps; {
    'E2_weights', 'baseline', {}
    'E2_weights', 'R = 0.001', {'R', 0.001}
    'E2_weights', 'R = 0.1', {'R', 0.1}
    'E2_weights', 'R = 1', {'R', 1}
    'E2_weights', 'Qf x 0.01', {'Qf', 0.01*Qf}
    'E2_weights', 'Qf x 0.1', {'Qf', 0.1*Qf}
    'E2_weights', 'Qf x 100', {'Qf', 100*Qf}
    'E2_weights', 'Q = 0', {'Q', zeros(4)}
    'E2_weights', 'Q theta heavy', {'Q', diag([5 0.1 50 1])}
    'E2_weights', 'old team weights', {'Q', diag([0.4 0.2 10 1]), 'R', 0.001, 'Qf', 1e5*diag([4 0.2 10 4])}}];
% 4.1 experiment 3: horizon length
for T = [1 1.5 2 2.5 3 4 6 8]
    exps(end+1,:) = {'E3_horizon', sprintf('T = %g s', T), {'T', T}};
end
% 4.1 experiment 4: initial guess
exps = [exps; {
    'E4_initial_guess', 'zeros', {}
    'E4_initial_guess', 'sine 10N 1Hz', {'U0', @(t) 10*sin(2*pi*t)}
    'E4_initial_guess', 'sine 30N 1Hz', {'U0', @(t) 30*sin(2*pi*t)}
    'E4_initial_guess', 'constant 40N', {'U0', @(t) 40 + 0*t}
    'E4_initial_guess', 'random 10N', {'U0', @(t) 10*randn(size(t))}
    'E4_initial_guess', 'zeros, theta0 = 150', {'x0', [0; 0; deg2rad(150); 0]}
    % 4.1 experiment 5: handling optimization failure (failing case, then change what the failure message points to)
    'E5_failure_handling', 'fail: Qf x 0.01', {'Qf', 0.01*Qf}
    'E5_failure_handling', 'fix: Qf x 0.1', {'Qf', 0.1*Qf}
    'E5_failure_handling', 'fail: T = 1.5', {'T', 1.5}
    'E5_failure_handling', 'fix: T = 1.5, q_p = 20', {'T', 1.5, 'Q', diag([20 0.1 1 0.1])}
    'E5_failure_handling', 'fail: Qf x 100', {'Qf', 100*Qf}
    'E5_failure_handling', 'try: Qf x 100, omega weight 1e3', {'Qf', diag([1e6 1e5 1e6 1e3])}
    'E5_failure_handling', 'fix: omega weight 10', {'Qf', diag([1e4 1e3 1e4 10])}
    'E5_failure_handling', 'fail: R = 0.001', {'R', 0.001}
    % 4.1 experiment 5: stuck at a low force limit, fixes that were tried
    'E5_low_force', 'Fmax = 20', {'Fmax', 20}
    'E5_low_force', 'Fmax = 20, T = 8', {'Fmax', 20, 'T', 8}
    'E5_low_force', 'Fmax = 20, sine guess', {'Fmax', 20, 'U0', @(t) 10*sin(2*pi*t)}
    'E5_low_force', 'Fmax = 20, R = 1', {'Fmax', 20, 'R', 1}
    'E5_low_force', 'Fmax = 20, R = 10', {'Fmax', 20, 'R', 10}}];

%% run everything
results = table();
for i = 1:size(exps,1)
    c = base;
    s = exps{i,3};
    for j = 1:2:numel(s)
        c.(s{j}) = s{j+1};
    end
    t = (0:round(c.T/Ts)-1)*Ts;
    rng(1);
    [X, U, K, cost, status] = ilqr(c.x0, c.U0(t), model, Ts, c.Q, c.R, c.Qf, c.Fmax);

    % check the nominal trajectory (same as in setup_cartpole.m)
    xN = X(:,end);
    xN(3) = wrap(xN(3));
    names = {'p_N', 'v_N', 'theta_N', 'omega_N'};
    vals = [xN(1) xN(2) rad2deg(xN(3)) xN(4)]; % angle in degrees for the table
    reason = {};
    for j = find(abs(xN) > cap)'
        reason{end+1} = sprintf('%s outside capture (%.2f)', names{j}, vals(j));
    end
    if max(abs(X(1,:))) > xmax
        reason{end+1} = sprintf('cart travel %.2f m', max(abs(X(1,:))));
    end
    ok = isempty(reason);
    reason = strjoin(reason, ', ');

    fprintf('%-34s it=%3d %-18s J=%-9.4g max|p|=%.2f  ok=%d  %s\n', exps{i,2}, ...
        numel(cost)-1, status, cost(end), max(abs(X(1,:))), ok, reason);

    % run it on the Simscape model
    Xnom = X; Unom = U; Knom = K; Fmax = c.Fmax; param.theta0 = rad2deg(c.x0(3));
    sim_out = sim('cartpole', 'StopTime', num2str(c.T + 4));
    xs = squeeze(sim_out.x.Data);
    balanced = norm([xs(end,1:2) wrap(xs(end,3)) xs(end,4)]) < 0.1;
    yesno = {'no', 'yes'};
    results = [results; table(string(exps{i,1}), string(exps{i,2}), string(yesno{ok+1}), string(reason), ...
        numel(cost)-1, string(status), round(cost(end)), round(xN(1),3), round(xN(2),3), round(rad2deg(xN(3)),1), ...
        round(xN(4),3), round(max(abs(X(1,:))),2), string(yesno{balanced+1}), ...
        'VariableNames', {'Experiment', 'Setting', 'Accepted', 'Reason', 'Iterations', 'Status', 'Cost', ...
        'p_end', 'v_end', 'theta_end_deg', 'omega_end', 'CartTravel_m', 'SimscapeBalanced'})];
    runs(i) = struct('X', X, 'U', U, 'cost', cost, 'ok', ok);
end
writetable(results, [out 'ilqr_experiment_summary.csv']);

% plots for the report: pole angle for the force and horizon experiments
for g = {'E1a_max_force', 'E3_horizon'}
    idx = find(strcmp(exps(:,1), g{1}));
    figure('Visible', 'off', 'Position', [100 100 800 450]); hold on;
    colors = [lines(7); 0 0 0]; % 8 different colors
    for n = 1:numel(idx)
        i = idx(n);
        style = '-'; if ~runs(i).ok, style = '--'; end % dashed = rejected
        plot((0:size(runs(i).X,2)-1)*Ts, rad2deg(unwrap(runs(i).X(3,:))), style, ...
            'Color', colors(n,:), 'LineWidth', 1.2, 'DisplayName', exps{i,2});
    end
    yline(0, 'k:', 'upright', 'HandleVisibility', 'off');
    grid on; xlabel('t [s]'); ylabel('\theta [deg]'); legend('Location', 'eastoutside');
    if strcmp(g{1}, 'E1a_max_force')
        title('Nominal pole angle for different force limits (dashed = rejected)');
    else
        title('Nominal pole angle for different horizons (dashed = rejected)');
    end
    saveas(gcf, [out 'figures/' g{1} '.png']);
end

%% Simscape: baseline with the right model, a wrong model, and no feedback
[Xnom, Unom, Knom] = ilqr(x0, zeros(1, 300), model, Ts, base.Q, base.R, base.Qf, base.Fmax);
Fmax = 40; param.theta0 = 180;
K_iLQR = Knom;
labels = {'tracking', 'tracking, 10% heavier pole + 2x cart friction', 'no feedback, same wrong model'};
for case_nr = 1:3
    param.m = 1; param.bc = 0.1; Knom = K_iLQR;
    if case_nr >= 2, param.m = 1.1; param.bc = 0.2; end % plant differs from the iLQR model
    if case_nr == 3, Knom = 0*K_iLQR; end
    sim_out = sim('cartpole', 'StopTime', '8');
    xs = squeeze(sim_out.x.Data);
    fprintf('%s: balanced at the end = %d\n', labels{case_nr}, norm([xs(end,1:2) wrap(xs(end,3)) xs(end,4)]) < 0.1);
end
param.m = 1; param.bc = 0.1; Knom = K_iLQR;

%% task 5.1: push on the cart, the supervisor in the controller block
% detects the fall, lets the pole hang still and swings it up again
fprintf('\npush on the cart (0.1 s at t = 6 s):\n');
figure('Visible', 'off', 'Position', [100 100 900 450]); hold on;
for dF_amp = [10 20 30 40 60 80]
    sim_out = sim('cartpole', 'StopTime', '20');
    xs = squeeze(sim_out.x.Data);
    tt = sim_out.x.Time;
    fell = any(abs(wrap(xs(tt > 6, 3))) > deg2rad(30));
    fprintf('%2d N: pole fell = %d, balanced at the end = %d\n', dF_amp, fell, ...
        norm([xs(end,1:2) wrap(xs(end,3)) xs(end,4)]) < 0.1);
    if dF_amp == 30 || dF_amp == 60 % one push the LQR handles, one that knocks the pole over
        plot(tt, rad2deg(abs(wrap(xs(:,3)))), 'LineWidth', 1.2, 'DisplayName', sprintf('%d N push', dF_amp));
    end
end
dF_amp = 0;
xline(6, 'r:', 'push', 'HandleVisibility', 'off');
grid on; xlabel('t [s]'); ylabel('|\theta| [deg] (0 = upright, 180 = down)'); ylim([0 190]); legend;
title('Push on the cart at t = 6 s (Simscape)');
saveas(gcf, [out 'figures/T51_disturbance_supervisor.png']);
