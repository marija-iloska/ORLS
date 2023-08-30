function [y, H, theta] = generate_data(T, dy, r,rt,  p_s, var_y)


% Choose random indices to be 0s
j = datasample(1:dy, round(p_s*dy));

% Generate random theta in the range between -rt, rt
theta = unifrnd(-rt,rt, dy, 1);

% Set chosen indices to 0s
theta(j) = 0;

% Create basis functions and data
H = unifrnd(-r, r, T, dy);
y = H*theta + mvnrnd(zeros(T,1), var_y*eye(T))';


end