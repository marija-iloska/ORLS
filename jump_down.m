<<<<<<< Updated upstream
function [theta_k, H, J,  Dk, k] = jump_down(y, k, Dk, theta_k, J, H, t, var_y)
=======
function [theta_k, H, J,  Dk, k, start] = jump_down(y, k, Dk, theta_k, J, H, t, t0, var_y, start)
>>>>>>> Stashed changes


for j = 1:k

    % Update current theta by jth basis function
<<<<<<< Updated upstream
    [theta_temp, D_temp, Hk_temp, J_temp, Hnew] = ols_downdates(y, theta_k, Dk, j, H, t, J);


    % COMPUTE time estimate    t ---> t+
    %[theta_temp, Sigma, J_temp] = time_update(y, Hk_temp, t, theta_temp, var_y, D_temp, J_temp);

=======
    [theta_temp, D_temp, Hk_temp,  Hnew] = ols_downdates(theta_k, Dk, j, H, t, start);

    % Compute predictive J(k,t) ---> J(k-1,t)
    [G, E, start_new] = pred_error_down(y, Hk_temp, H, t, t0, var_y, J, start, j, D_temp, theta_temp);
    J_temp = J - (G*G' + 2*G*E);
>>>>>>> Stashed changes

    % Corresponding variables store
    H_store{j} = Hnew;
    D_store{j} = D_temp;
    theta_store{j} = theta_temp;
    J_store(j) = J_temp;
<<<<<<< Updated upstream
    %Sigma_store{j} = Sigma;
=======
    start_store{j} = start_new;
>>>>>>> Stashed changes
end

% Choose min J to update
min_idx = find(J_store == min(J_store));
% Ws = exp(-(J_store-min(J_store)));
% Ws = Ws./sum(Ws);
% min_idx = datasample(1:k, 1, 'Weights', Ws);

% Update all parameters
theta_k = theta_store{min_idx};
H  = H_store{min_idx};  % For optimization I could save indices here
J = J_store(min_idx);
Dk = D_store{min_idx};
<<<<<<< Updated upstream
%Sigma = Sigma_store{min_idx};
%Dk = Sigma/var_y;
=======
start = start_store{min_idx};

>>>>>>> Stashed changes
k = length(theta_k);

end