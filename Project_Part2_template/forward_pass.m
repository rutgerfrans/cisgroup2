function [Xnew,Unew] = forward_pass(x0,X,U,d,K,alpha,model,Ts,Fmax)
% X, U: nominal trajectory
% d: feedforward corrections (1,N)--> this corresponds to k in the class
% notes. Here denoted as d to avoid conflict with time index. 
% K: feedback gains of size (1,4,N)
% alpha: line-search step size  unew= unom +\alpha*k - K (xnew-xnom)

N = numel(U);

% Initialization
Xnew = zeros(4,N+1);
Unew = zeros(1,N);
Xnew(:,1) = x0;


%%%% Implement forward pass here %%%

for k = 1:N
    Unew(:, k) = U(:, k) + alpha*d(:, k) - K(:, :, k)*(Xnew(:,k) - X(:, k));
    Unew(:, k) = max(-Fmax,min(Fmax, Unew(:, k)));
    Xnew(:, k+1) = discrete_step(Xnew(:, k), Unew(:, k), model, Ts);
end

end