function [G, E, start] = pred_error_down(y, Hk, H,  t, t0, var_y, J_old, start, j, D, theta)

% Get k
k = length(Hk(1,1:end-1));


 % k x k at t0
Dk_test = inv(Hk(1:t0, 1:k)'*Hk(1:t0, 1:k));
theta_k_test = Dk_test*Hk(1:t0, 1:k)'*y(1:t0);

[theta_kk, Dkk] = start{:};


% k+1  x  k+1 at t0
[theta_k, Dk, ~, ~] = ols_downdates(theta_kk, Dkk, j, H, t);

% Starts
start = {theta_k, Dk};


% Initialize
G = [];
THETA = [];

for i = t0+1:t

    % Compute Ai
    A = Hk(i, 1:k)*Dk*Hk(1:i-1, 1:k)'*Hk(1:i-1, k+1) - Hk(i,k+1);

    % Compute Gi
    G(end+1) = A*theta_kk(end);

    % Stack theta estimates into matrix
    THETA = [THETA; theta_k'];

    if (i == t)
        break
    end


    % Compute theta_(k+1, t-1), check Dk indices
    [theta_kk, Dkk, ~] = time_update(y, Hk(1:i, :), i, theta_kk, var_y, Dkk, J_old);
    [theta_k, Dk, ~] = time_update(y, Hk(1:i, 1:k), i, theta_k, var_y, Dk, J_old);

end

% Compute predictive residual error
E = y(t0+1:t) - sum( Hk(t0+1:t, 1:k).*THETA , 2 );





end