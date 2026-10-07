function J = trajectory_cost(X,U,Q,R,Qf)
% J is the total cost corresponding to the control trajectory U from
% initial state x0 
% Use the objective function given in the project manual


N = numel(U);

%%%% Implement cost calculation here %%%%

J = 1/2 * X(:, N + 1)' * Qf * X(:, N+1);
sumCost = 0;
for k = 1:N
    sumCost = sumCost + X(:, k)' * Q * X(:, k) + U(:, k)' * R * U(:, k);
end
J = J + 1/2 * sumCost;

end