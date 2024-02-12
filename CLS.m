clc; clear; close all; 

t = [0:0.1:10]'; % time vector 10 seconds
data = load("data4.mat"); % load data
y_m = data.y_m; % measurements

% Case 1
y_m_1 = y_m(1:100);
y_m_2 = y_m(101);
% Basis Functions
H1 = [t(1:100).^2 sin(t(1:100)) cos(t(1:100)) exp(t(1:100))]; % basis function matrix
H2 = [t(101).^2 sin(t(101)) cos(t(101)) exp(t(101))]; % basis function matrix
% Get Constrained Estimate
xe_bar_1 = ((H1'*H1)^(-1)*H1'*y_m_1)';
K = (H1'*H1)^(-1)*H2'*(H2*(H1'*H1)^(-1)*H2')^(-1); % gain
xe_1 = (xe_bar_1'+K*(y_m_2-H2*xe_bar_1'))
y_e_1 = xe_1(1)*t.^2 + xe_1(2)*sin(t) + xe_1(3)*cos(t) + xe_1(4)*exp(t); % estimation
norm_e_1 = norm(y_m(97:101) - y_e_1(97:101)); % residual error

% Case 2
y_m_1 = y_m(1:99);
y_m_2 = y_m(100:101);
% Basis Functions
H1 = [t(1:99).^2 sin(t(1:99)) cos(t(1:99)) exp(t(1:99))]; % basis function matrix
H2 = [t(100:101).^2 sin(t(100:101)) cos(t(100:101)) exp(t(100:101))]; % basis function matrix
% Get Constrained Estimate
xe_bar_2 = ((H1'*H1)^(-1)*H1'*y_m_1)';
K = (H1'*H1)^(-1)*H2'*(H2*(H1'*H1)^(-1)*H2')^(-1); % gain
xe_2 = (xe_bar_2'+K*(y_m_2-H2*xe_bar_2'))
y_e_2 = xe_2(1)*t.^2 + xe_2(2)*sin(t) + xe_2(3)*cos(t) + xe_2(4)*exp(t); % estimation
norm_e_2 = norm(y_m(97:101) - y_e_2(97:101)); % residual error

% Case 3
y_m_1 = y_m(1:98);
y_m_2 = y_m(99:101);
% Basis Functions
H1 = [t(1:98).^2 sin(t(1:98)) cos(t(1:98)) exp(t(1:98))]; % basis function matrix
H2 = [t(99:101).^2 sin(t(99:101)) cos(t(99:101)) exp(t(99:101))]; % basis function matrix
% Get Constrained Estimate
xe_bar_3 = ((H1'*H1)^(-1)*H1'*y_m_1)';
K = (H1'*H1)^(-1)*H2'*(H2*(H1'*H1)^(-1)*H2')^(-1); % gain
xe_3 = (xe_bar_3'+K*(y_m_2-H2*xe_bar_3'))
y_e_3 = xe_3(1)*t.^2 + xe_3(2)*sin(t) + xe_3(3)*cos(t) + xe_3(4)*exp(t); % estimation
norm_e_3 = norm(y_m(97:101) - y_e_3(97:101)); % residual error

% Case 4
y_m_1 = y_m(1:97);
y_m_2 = y_m(98:101);
% Basis Functions
H1 = [t(1:97).^2 sin(t(1:97)) cos(t(1:97)) exp(t(1:97))]; % basis function matrix
H2 = [t(98:101).^2 sin(t(98:101)) cos(t(98:101)) exp(t(98:101))]; % basis function matrix
% Get Constrained Estimate
xe_bar_4 = ((H1'*H1)^(-1)*H1'*y_m_1)';
K = (H1'*H1)^(-1)*H2'*(H2*(H1'*H1)^(-1)*H2')^(-1); % gain
xe_4 = (xe_bar_4'+K*(y_m_2-H2*xe_bar_4'))
y_e_4 = xe_4(1)*t.^2 + xe_4(2)*sin(t) + xe_4(3)*cos(t) + xe_4(4)*exp(t); % estimation
norm_e_4 = norm(y_m(97:101) - y_e_4(97:101)); % residual error

% Case 4 check
K_check = H2^(-1);
xe_4_check = K_check * y_m_2;

% plot residual error norm history
norms = [norm_e_1 norm_e_2 norm_e_3 norm_e_4];
cases = ["case 1" "case 2" "case 3" "case 4"];

figure;
bar(cases,norms)
xlabel('Case')
ylabel('norm residual error')