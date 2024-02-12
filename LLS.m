clc; clear; close all; 

t = [0:0.1:10]'; % time vector 10 seconds
data = load("data4.mat"); % load data
y_m = data.y_m; % measurements

% model 1
H_1 = [t sin(t) cos(t) exp(t)]; % basis function matrix

a = (H_1'*H_1)^(-1)*H_1'*y_m % find coefficients
y_e_1 = a(1)*t + a(2)*sin(t) + a(3)*cos(t) + a(4)*exp(t); % estimation
e_1 = y_m - y_e_1; % residual error
mu_e_1 = mean(e_1); % mean of residual error
sig_e_1 = std(e_1); % standard deviation of residual error
fprintf('model 1 mean residual error: %f', mu_e_1)
fprintf('\nmodel 1 standard deviation residual error: %f', sig_e_1)

% plot estimation with measurement
figure;
plot(t,y_m,t,y_e_1)
xlabel('time')
ylabel('y')
legend('Measurement','Estimate');
% plot residual error
figure;
plot(t,e_1,'*')
xlabel('time')
ylabel('residual error')
%model 2
H_2 = [t.^2 sin(t) cos(0.9*t) exp(t)]; % basis function matrix

b = (H_2'*H_2)^(-1)*H_2'*y_m % find coefficients
y_e_2 = b(1)*t.^2 + b(2)*sin(t) + b(3)*cos(0.9*t) + b(4)*exp(t); % estimation
e_2 = y_m - y_e_2; % residual error
mu_e_2 = mean(e_2); % mean of residual error
sig_e_2 = std(e_2); % standard deviation of residual error
fprintf('\nmodel 2 mean residual error: %f', mu_e_2)
fprintf('\nmodel 2 standard deviation residual error: %f', sig_e_2)
% plot estimation with measurement
figure;
plot(t,y_m,t,y_e_2)
xlabel('time')
ylabel('y')
legend('Measurement','Estimate');
% plot residual error
figure;
plot(t,e_2,'*')
xlabel('time')
ylabel('residual error')
% model 3
H_3 = [t.^2 sin(t) cos(t) exp(0.9*t)]; % basis function matrix

c = (H_3'*H_3)^(-1)*H_3'*y_m % find coefficients
y_e_3 = c(1)*t.^2 + c(2)*sin(t) + c(3)*cos(t) + c(4)*exp(0.9*t); % estimation
e_3 = y_m - y_e_3; % residual error
mu_e_3 = mean(e_3); % mean of residual error
sig_e_3 = std(e_3); % standard deviation of residual error
fprintf('\nmodel 3 mean residual error: %f', mu_e_3)
fprintf('\nmodel 3 standard deviation residual error: %f', sig_e_3)
% plot estimation with measurement
figure;
plot(t,y_m,t,y_e_3)
xlabel('time')
ylabel('y')
legend('Measurement','Estimate');
% plot residual error
figure;
plot(t,e_3,'*')
xlabel('time')
ylabel('residual error')
% model 4
H_4 = [t.^2 sin(t) cos(t) exp(t)]; % basis function matrix

d = (H_4'*H_4)^(-1)*H_4'*y_m % find coefficients
y_e_4 = d(1)*t.^2 + d(2)*sin(t) + d(3)*cos(t) + d(4)*exp(t); % estimation
e_4 = y_m - y_e_4; % residual error
mu_e_4 = mean(e_4); % mean of residual error
sig_e_4 = std(e_4); % standard deviation of residual error
fprintf('\nmodel 4 mean residual error: %f', mu_e_4)
fprintf('\nmodel 4 standard deviation residual error: %f\n', sig_e_4)
% plot estimation with measurement
figure;
plot(t,y_m,t,y_e_4)
xlabel('time')
ylabel('y')
legend('Measurement','Estimate');
% plot residual error
figure;
plot(t,e_4,'*')
xlabel('time')
ylabel('residual error')