function [X,U,K,costHistory] = ilqr(x0,U0,model,Ts,Q,R,Qf,Fmax)

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
error_threshold = 0.001;

costHistory = [trajectory_cost(X, U, Q, R, Qf)];
cost_diff = 10^9;
iter = 1;
while cost_diff > error_threshold && iter < 100
    [A, B] = linearize_trajectory(X, U, model, Ts);
    [d, K] = backward_lqr(X, U, A, B, Q, R, Qf);

    not_improved = true;
    alpha = 1;
    iter = 1;
    while not_improved && iter < 100
        [Xnew, Unew] = forward_pass(x0,X,U,d,K,alpha,model,Ts,Fmax);
        Jnew = trajectory_cost(Xnew, Unew, Q, R, Qf);
        if Jnew < costHistory(:, end)
            not_improved = false;
            cost_diff = costHistory(:, end) - Jnew;
            costHistory = [costHistory, Jnew];
            X = Xnew;
            U = Unew;
        else
            alpha = beta * alpha;
            iter = iter + 1;
        end
    end
end
[d, K] = backward_lqr(X, U, A, B, Q, R, Qf);

end