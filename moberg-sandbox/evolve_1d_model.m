% modified wilson-cowan flicker phosphene model
clear; clc; 
rng(1);

% spatial initialization
Ldom = 50; % length of domain 
N = 128; % number of points
dx = Ldom/N; % spacing
x = (0:N-1)*dx - Ldom/2; % x vector, shift over so it is centered at zero

% time (ms)
dt = 0.1; % timestep spacing
tEnd = 2000; % final time
t = 0:dt:tEnd; % time vector
M = numel(t); % total number of time points

% inline functions
F = @(u) 1./(1+exp(-u)); % sigmoid F
Kern = @(x,s) exp(-x.^2/s^2)/(s*sqrt(pi)); % kernel
H = @(y) double(y>0); % step function

% parameters
tau_e=10; 
tau_i=20;
a_ee=10; 
a_ei=12; 
a_ie=8.5; 
a_ii=3;
th_e=2; % theta_e 
th_i=3.5; % theta_i
sig_e=2; 
sig_i=2.5*sig_e;
A=0.8; 
T=108; % forcing function period                  
g_e=A; 
g_i=0;                   

% kernels - same each time so only need to compute once
% use ifftshift to move the peak of the kernel from index N/(2+1) to index 1
% compute fft to transform function into vector
Ke_hat = fft(ifftshift(Kern(x,sig_e))); 
Ki_hat = fft(ifftshift(Kern(x,sig_i)));
% verify kernel normalization and zero-frequency response
sum(Kern(x,sig_e))*dx % ≈ 1: kernel normalized
real(Ke_hat(1))*dx  % ≈ 1: uniform fields pass through unchanged
max(abs(imag(Ke_hat)))  % ≈ 0: confirms correct centering
       

% initial condition: near-uniform plus small noise
Ue = 0.1 + 0.01*randn(1,N);
Ui = 0.1 + 0.01*randn(1,N);

% store Ue and Ui
hist_e = zeros(M,N);
hist_i = zeros(M,N);

% use fft to transform Ue and Ui
% shifted kernel real and even so fft is real
% apply real() to remove roundoff errors
% compute Ue and Ui using FW Euler
for k = 1:M
    S = H(sin(2*pi*t(k)/T) - 0.8);   % A is applied through g_e below
    KeUe = dx*real(ifft(fft(Ue).*Ke_hat)); 
    KiUi = dx*real(ifft(fft(Ui).*Ki_hat));
    Ue = Ue + dt/tau_e*(-Ue + F(a_ee*KeUe - a_ie*KiUi - th_e + g_e*S));
    Ui = Ui + dt/tau_i*(-Ui + F(a_ei*KeUe - a_ii*KiUi - th_i + g_i*S));
    hist_e(k,:) = Ue;
    hist_i(k,:) = Ui;
end

figure;
imagesc(x, t(t>=1500), hist_e(t>=1500,:)); axis xy; colorbar
xlabel('x'); ylabel('t (ms)'); title(sprintf('U_e, T = %g ms', T))

figure;
imagesc(x, t(t>=1500), hist_i(t>=1500,:)); axis xy; colorbar
xlabel('x'); ylabel('t (ms)'); title(sprintf('U_i, T = %g ms', T))