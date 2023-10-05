clear all
close all
clc

% Settings
var_y = 1;     % Variance
ps = 20;       % Sparsity percent
dy = 50;      % System dimension
T = 100;      % Time series length
r = 1;         % Range of input data H
rt = 3;      % Range of theta
t0 = 50;

%Create data
[y, H, theta] = generate_data(T, dy, r, rt,  ps, var_y);      
idx_h = find(theta ~= 0)';

% Define initial batch
y0 = y(1:t0);
H0 = H(1:t0, :);

% Define initial batch terms
xy0 = H0'*y0;
xx0 = H0'*H0;

% EIG
a = eig(xx0);
step = 0.01*t0/max(real(a));

% Initial estimate
[B, STATS] = lasso(H0, y0, 'CV', 5);
theta_olin = B(:, STATS.IndexMinMSE);
epsilon = 1e-3;

% Initialize terms
xy = zeros(dy,1);
xx = zeros(dy,dy);

J_pred = [];
J_incr = 0;


for t = t0+1:T-1

    % Updates
    xx = xx + H(t,:)'*H(t,:);
    xy = xy + H(t,:)'*y(t);
    
    [theta_olin, loss{t}] = olin_lasso(xy0, xx0, xy, xx, theta_olin, epsilon, step, t0, t, dy);

    J_pred(end+1) = sum( (y(1:t+1) - H(1:t+1,:)*theta_olin).^2);
    J_incr(end+1) = J_incr(end) + (y(t+1) - H(t+1, :)*theta_olin)^2;
end


plot(J_incr)
hold on
plot(J_pred)


