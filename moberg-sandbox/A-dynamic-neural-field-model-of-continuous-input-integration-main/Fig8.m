%% This code reproduces Fig.8 (a) from the article 
%  'A dynamic neural field model of continuous input integration'
%  by W. Wojtak et al.

clear; clc
% spatial stepping
L = 4*pi; % domain length
N = 2^12; % number of spatial points
dx = 2*L/N; % spatial step
xDim = (-L+(0:N-1)*dx); % spatial vector

% timestepping
dt = .01; 
tspan = 0:dt:50; % create time vector to 50 timesteps
M = numel(tspan); % length of time vector

% define inline functions
kernel = @(x,A_ex,s_ex,A_in,s_in,g_i) A_ex* exp(-0.5 * (x).^2 / s_ex^2) - A_in * exp(-0.5 * (x).^2 / s_in^2) - g_i; % mexican hat kernel
gauss = @(x,mu,sigma) exp(-0.5*(x-mu).^2/sigma^2); % generic gaussian
sigmoid = @(x,beta,theta)  1./ (1 + exp(-beta*(x-theta))); % sigmoid

% params
p(1) = 1000;   % beta
p(2) = 2;      % A_ex   
p(3) = 1.25;   % s_ex   
p(4) = 1;      % A_inh
p(5) = 2.5;    % s_inh
p(6) = 0.1;    % w_inh
p(7) = 0.5;    % theta

beta = p(1); theta = p(7); tau = 1;

% initial conditions
K = 0.0;
u_field =  -theta * ones(1, N); % u initial condition
v_field = K - u_field; % v initial condition
history_u = zeros(M, N); % u history array
history_v = zeros(M, N); % v history arary

% connectivity function
w = kernel(xDim,p(2),p(3),p(4),p(5),p(6)); % construct mexican hat kernel using specified 
wHat = fft(w);

%% Input
A_I = 1; 
sigma_I = 1;
Input = zeros(M, N);
Input_pattern = A_I * gauss(xDim, 0, sigma_I);
Input(1/dt:2/dt, :) = repmat(Input_pattern,1+(1/dt),1);

%% Main loop
for i = 1:M
f = sigmoid(u_field, beta, theta);
convolution = dx * ifftshift(real(ifft(fft(f) .* wHat)));
% 2-field model: time integration with FW Euler
u_field = u_field + dt/tau * (-u_field + convolution + v_field + Input(i, :)); 
v_field = v_field + dt/tau * (-v_field - convolution + u_field);
history_u(i,:) = u_field; history_v(i,:) = v_field;
end


%% Plot results
figure
plot(xDim,u_field,'k','linewidth',3), hold on
plot(xDim,v_field,'--k','linewidth',3), hold on
plot(xDim,theta*ones(1,N),':k','linewidth',2);
xlabel('x'); ylabel('u(x), v(x)');
ax = gca;
set(gca,  'FontSize', 20)
set(gca,  'XLim', [-L L])

