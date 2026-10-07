function [X,U,K,costHistory,status] = ilqr(x0,U0,model,Ts,Q,R,Qf,Fmax)

x0 = x0(:);
U0 = min(Fmax,max(-Fmax,U0(:, :)));
N = numel(U0);

% Initial trajectory
[X,U] = forward_pass(x0,zeros(4,N+1),U0,zeros(1,N),zeros(1,4,N), 0,model,Ts,Fmax);

%%% Implement iLQR here %%%
%%%%% Make sure the implementation has the following:
%%%%% Line search (iterate over different values of alpha in the forward pass)
%%%%% Stopping and convergence criteria

beta = 0.5;
rel_threshold = 1e-4; % stop when cost decreases less than 0.01%
max_iter = 300;
max_ls = 30;

costHistory = trajectory_cost(X, U, Q, R, Qf);
status = 'max iterations';
for iter = 1:max_iter
    [A, B] = linearize_trajectory(X, U, model, Ts);
    [d, K] = backward_lqr(X, U, A, B, Q, R, Qf);

    % line search, halve alpha until the cost goes down
    alpha = 1;
    improved = false;
    for ls = 1:max_ls
        [Xnew, Unew] = forward_pass(x0,X,U,d,K,alpha,model,Ts,Fmax);
        Jnew = trajectory_cost(Xnew, Unew, Q, R, Qf);
        if Jnew < costHistory(end)
            improved = true;
            break
        end
        alpha = beta * alpha;
    end
    if ~improved
        status = 'line search failed';
        break
    end

    cost_diff = costHistory(end) - Jnew;
    costHistory = [costHistory, Jnew];
    X = Xnew;
    U = Unew;

    if cost_diff < rel_threshold * Jnew
        % a tiny accepted step means we are stuck (the backwards pass
        % barely contributes to the upgrade)
        if alpha >= 1/16
            status = 'converged';
        else
            status = 'stalled';
        end
        break
    end
end

% final gains, linearized along the accepted trajectory
[A, B] = linearize_trajectory(X, U, model, Ts);
[~, K] = backward_lqr(X, U, A, B, Q, R, Qf);

end
