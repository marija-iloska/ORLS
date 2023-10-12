function [theta_k, Dk, Hk, J, H] = ols_updates(y, H, k, j, t, t0, var_y, Dk, theta_k, J_old)

                                 
    % Current input data
    K = length(H(1,:));
    Hk = H(1:(t-1), 1:k);
    Dk_old = Dk;
    y_past = y(1:(t-1));

    Pk_norm = eye(t-1) - Hk*Dk*Hk';

    % Take the new observation column h(k+1)
    hk = H(1:(t-1),k+j);

    % Reuseable terms
    temp = hk'*Pk_norm;
    xd = temp*y_past;
    d = Dk*Hk'*hk;

    % Compute terms of D(k+1)
    DK22 = 1 / ( temp*hk );
    DK12 = - d*DK22;
    DK21 = DK12';
    DK11 = Dk + DK22*(d*d');

    % Update D(k) to D(k+1)
    Dk = [ DK11, DK12 ;  DK21, DK22 ];
    Dkk_old = Dk;

    % Update theta_k
    theta_k = [ theta_k - d*xd*DK22;    xd*DK22 ];

    % Update Hk in time and k
    Hk = H(1:t, [1:k, k+j]);

    % Update original available data H (swap column order )
    H = H(:, [1:k, k+j, setdiff( (k+1):K, (k+j) ) ]);


    % Compute Jk ---> Jk+
    %J_old = J_old - theta_k(end)^2/DK22;


    % Compute Jk ---> Jk+
    %J =  (y(t) - Hk(t, :)*theta_k)^2;
    %J = sum( (y(1:t) - Hk*theta_k).^2);
    [G, V] = pred_error(y, Hk, t, t0, var_y, J_old);
    J = J_old + (G*G' + 2*G*V);
    %J =  G*G' + 2*G*V ;



end
