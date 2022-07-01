function [theta_k, Dk, Hk] = ols_downdates(theta_k, Dk, min_k, H, t)


% Code for D change (before update)
% Let be the column removed (or min_k)
K = length(Dk(1,:));
idx = 1:K;
idx = setdiff(idx, min_k);


% Hk update
Hk = H(1:t, idx);

Dknew = Dk(idx, idx);
Dknew(K, 1:K-1) = Dk(min_k, idx);
Dknew(:, K) = Dk([idx, min_k], min_k);

% Now Dk downdate code
DK11 = Dknew(1:K-1, 1:K-1);
DK12 = Dknew(1:K-1, K);
%DK21 = Dknew(K, 1:K-1);
DK22 = Dknew(K,K);

Dk = DK11 - DK12*DK12' /( DK22^2 );


idx = 1:length(theta_k);
idx = setdiff(idx, min_k);

% Ratio
ratio = DK12/DK22;

% Update rest
theta_k(idx) = theta_k(idx) - ratio*theta_k(min_k);

% Rest
theta_k(min_k) = [];







end