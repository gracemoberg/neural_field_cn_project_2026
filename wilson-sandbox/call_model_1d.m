% modified wilson-cowan flicker phosphene model on 1d spatial domain
clear; clc;

% access default parameters unless otherwise specified
par = default_par();
par.tEnd = 6000; % final time
t_plot = 5000; % set portion of run to plot

% forcing parameters
A = 0.8; % forcing amplitude
T = 42; % forcing function period
q = 0; % feed-forward inhibition fraction (g_i = q*A in the paper)

H = @(y) double(y>0); % step function

% can define different forcing functions
% each returns the drive at time t (scalar or 1xN)
S0 = @(t) A*H(sin(2*pi*t/T) - 0.8); % paper's S(t) = A*H(sin(2*pi*t/T)-0.8)

% forcing handles take (t,x) so spatial forcing is possible
forcing_e = @(t,x) S0(t); % drive to E population
forcing_i = @(t,x) q*S0(t); % drive to I population
% spatial forcing example:
% forcing_e = @(t,x) S0(t).*(1 + 0.1*cos(2*pi*4*x/par.Ldom));

% grid (same construction as inside run_flicker, needed here to build ICs)
dx = par.Ldom/par.N;
x = (0:par.N-1)*dx - par.Ldom/2;

% initial condition: near-uniform plus small noise
rng(7);
% Ue0 = 0.1 + 0.01*randn(1,par.N);
% Ui0 = 0.1 + 0.01*randn(1,par.N);
% other options:
Ue0 = rand(1,par.N);  
Ui0 = rand(1,par.N); % large random
 %Ue0 = 0.1*ones(1,par.N) + 0.05*exp(-x.^2/4); % localized bump on uniform background
%Ue0 = 0.1 + 0.01*sin(2*pi*x/par.Ldom); % single seeded wavenumber (k = 4)
 %Ui0 = Ue0; % same state in I as in E

[x,t,hist_e,hist_i] = run_flicker_1d(forcing_e, forcing_i, par, Ue0, Ui0);

figure;
imagesc(x, t(t>=t_plot), hist_e(t>=t_plot,:)); axis xy; colorbar
xlabel('x'); ylabel('t (ms)'); title(sprintf('U_e, T = %g ms', T))
figure;
imagesc(x, t(t>=t_plot), hist_i(t>=t_plot,:)); axis xy; colorbar
xlabel('x'); ylabel('t (ms)'); title(sprintf('U_i, T = %g ms', T))

% pattern measure D from the paper over the last 10 periods (D = 0 means no pattern)
dx = x(2) - x(1);
last = t > par.tEnd - 10*T;
D = mean(sum(abs(hist_e(last,:) - hist_e(last,1)), 2)) * dx % average over last ten periods, diff between all spatial points and 0 spatial point 

