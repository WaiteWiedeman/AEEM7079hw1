clc; clear; close all; 

t = [0:0.1:10]'; % time vector 10 seconds
data = load("data4.mat"); % load data
y_m = data.y_m; % measurements

W = ones(1,98);
w_i = [1 1e2 1e4 1e6 1e8];

% weight 1
W_1 = diag([W w_i(1) w_i(1) w_i(1)]);
[d_1, y_e_1] = estimate(W_1,y_m,t);
norm_e_1 = norm(y_m(98:101) - y_e_1(98:101)); % residual error

fprintf('weight 1 norm residual error: %f\n', norm_e_1)
fprintf('weight 1 estimation coefficients: %f, %f, %f, %f\n', d_1)

% weight 2
W_2 = diag([W w_i(2) w_i(2) w_i(2)]);
[d_2, y_e_2] = estimate(W_2,y_m,t);
norm_e_2 = norm(y_m(98:101) - y_e_2(98:101)); % residual error

fprintf('weight 2 norm residual error: %f\n', norm_e_2)
fprintf('weight 2 estimation coefficients: %f, %f, %f, %f\n', d_2)

% weight 3
W_3 = diag([W w_i(3) w_i(3) w_i(3)]);
[d_3, y_e_3] = estimate(W_3,y_m,t);
norm_e_3 = norm(y_m(98:101) - y_e_3(98:101)); % residual error

fprintf('weight 3 norm residual error: %f\n', norm_e_3)
fprintf('weight 3 estimation coefficients: %f, %f, %f, %f\n', d_3)

% weight 4
W_4 = diag([W w_i(4) w_i(4) w_i(4)]);
[d_4, y_e_4] = estimate(W_4,y_m,t);
norm_e_4 = norm(y_m(98:101) - y_e_4(98:101)); % residual error

fprintf('weight 4 norm residual error: %f\n', norm_e_4)
fprintf('weight 4 estimation coefficients: %f, %f, %f, %f\n', d_4)

% weight 5
W_5 = diag([W w_i(5) w_i(5) w_i(5)]);
[d_5, y_e_5] = estimate(W_5,y_m,t);
norm_e_5 = norm(y_m(98:101) - y_e_5(98:101)); % residual error

fprintf('weight 5 norm residual error: %f\n', norm_e_5)
fprintf('weight 5 estimation coefficients: %f, %f, %f, %f\n', d_5)

% plot residual error norm history
norms = [norm_e_1 norm_e_2 norm_e_3 norm_e_4 norm_e_5];
weights = ["w_1" "w_2" "w_3" "w_4" "w_5"];

figure;
bar(weights,norms)
xlabel('weight')
ylabel('norm residual error')

function [d, y_e] = estimate(W,y_m,t)
H = [t.^2 sin(t) cos(t) exp(t)]; % basis function matrix
d = (H'*W*H)^(-1)*H'*W*y_m; % find coefficients
y_e = d(1)*t.^2 + d(2)*sin(t) + d(3)*cos(t) + d(4)*exp(t); % estimation
end