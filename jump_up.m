<<<<<<< Updated upstream
function [theta_k, H, J,  Dk, k] = jump_up(y, dx, k, Dk, theta_k, J, H, t)
=======
function [theta_k, H, J,  Dk, k, start] = jump_up(y, dx, k, Dk, theta_k, J, H, t, t0, var_y, start)

>>>>>>> Stashed changes


for j = 1:(dx - k)

    % Update current theta by jth basis function
<<<<<<< Updated upstream
    [theta_temp, D_temp, ~, J_temp, Hnew] = ols_updates(y, H, k, j, t, Dk, theta_k, J);
    if(J_temp < 0)
        J_temp
    end
=======
    [theta_temp, D_temp, Hk_temp,  Hnew] = ols_updates(y, H, k, j, t, Dk, theta_k);

    % Compute J(k,t) ---> J(k+1, t)
    [G, E, start_new] = pred_error_up(y, Hk_temp, t, t0, var_y, J, start, D_temp, theta_temp);
    %J_temp = J - (G*G' - 2*G*E);
    J_temp = J + (G*G' + 2*G*E);

>>>>>>> Stashed changes
                           
    % Corresponding variables store
    H_store{j} = Hnew;
    D_store{j} = D_temp;
    theta_store{j} = theta_temp;
    J_store(j) = J_temp;
    start_store{j} = start_new;


end

% Choose min J to update
min_idx = find(J_store == min(J_store));
% Ws = exp(-(J_store-min(J_store)));
% Ws = Ws./sum(Ws);
% min_idx = datasample(1:dx-k, 1, 'Weights', Ws);

% Update all parameters
if (isempty(min_idx))
    disp('stop')
end
theta_k = theta_store{min_idx};
H  = H_store{min_idx};  % For optimization I could save indices here
J = J_store(min_idx);
Dk = D_store{min_idx};
k = length(theta_k);
start = start_store{min_idx};

end