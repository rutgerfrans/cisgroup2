function captured = test_capture(x0, K_lqr, model, Ts, Fmax, xmax, Tsim)
% returns variable *captured as true if the LQR stabilizes the cart pole
% from the initial x0
% used to iterate over many starting positions to find pole capture zone

x = x0;
N = round(Tsim/Ts);
captured = false;
for k = 1 : N % repeat the simulink steps at each loop iteraiton
    u = - K_lqr * x; %  LQR state feedback formula
    u = min(max(u, -Fmax), Fmax); % force limit implementation
    x = discrete_step(x, u, model, Ts); % one sample step
    if abs(x(1)) > xmax || abs(x(3)) > pi/2 % cart limit hit || pole fell down
        return
    end
end
captured = norm(x) < 0.1; % ok if the pole settles 10 degrees near upright
end