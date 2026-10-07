function [d,K] = backward_lqr(X,U,A,B,Q,R,Qf)
% A is of size (4,4,N), A matrices corresponding to linearized trajectory
% B is of size (4,1,N), B matrices correspondng to linearized trajectory


N = numel(U);  %number of time steps

%%% Initialization
d = zeros(1,N); % feedforward terms 
% Note: d corresponds to the feedforward k term in the class notes

K = zeros(1,4,N); %state feedback gains

%%% Implement backward pass here 

S = Qf;
s = Qf * X(:, N + 1);

for i = N:-1:1
    Huu_inv = (R + B(:,:,i)'*S*B(:,:,i))^-1;
    K(:, :, i) = Huu_inv*B(:,:,i)'*S*A(:,:,i);
    d(:, i) = -Huu_inv*(R*U(:,i) + B(:,:,i)'*s);
    s = Q * X(:,i) + A(:,:,i)'*s + A(:,:,i)'*S*B(:,:,i)*d(:,i);
    S = Q + A(:,:,i)'*S*(A(:,:,i) - B(:,:,i)*K(:,:,i));
end

end