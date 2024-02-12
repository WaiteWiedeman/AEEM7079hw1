clc; clear; close all; 

t = [0:0.1:10]'; % time vector 10 seconds
data = load("data4.mat"); % load data
y_m = data.y_m; % measurements
m = length(y_m); % number of measurements

% Weight and H Matrix
std = 0.1;
W_k = std^-2;
H = [t.^2 sin(t) cos(t) exp(t)]; % basis function matrix

% Initial Conditions for Sequential Algorithm
alpha = 1e3;
beta=[1e-2 1e-3 1e-4 1e-4]';
P1 = inv(1/alpha/alpha*eye(4) + H(1,:)'*W_k*H(1,:));
x1 = P1*(1/alpha*beta + H(1,:)'*W_k*y_m(2));

% Sequential Least Squares
xk = zeros(m-1,4); xk(1,:) = x1';
p = zeros(m-1,4); p(1,:) = diag(P1)'; pp = P1;
for i=1:m-2
 k = pp*H(i+1,:)'*inv(H(i+1,:)*pp*H(i+1,:)'+inv(W_k));
 pp = (eye(4)-k*H(i+1,:))*pp;
 xk(i+1,:) = xk(i,:)+(k*(y_m(i+2)-H(i+1,:)*xk(i,:)'))';
 p(i+1,:) = diag(pp)';
end

% Plot Results
figure;
plot(t(1:100),xk(:,1),t(1:100),xk(:,2),t(1:100),xk(:,3),t(1:100),xk(:,4))
xlabel('time')
ylabel('x_k')
legend('x_k 1','x_k 2','x_k 3','x_k 4');

figure;
title('Diagonal Elements of P_k')
plot(t(1:100),p(:,1),t(1:100),p(:,2),t(1:100),p(:,3),t(1:100),p(:,4))
xlabel('time')
ylabel('P_k')
legend('P_k 1','P_k 2','P_k 3','P_k 4');