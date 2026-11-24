% deriving spatially homogeneous (stable) equilibrium in the unforced
% system

par = default_par(); % from shared files folder
dx = par.Ldom/par.N; % spatial domain 
x = (0:par.N-1)*dx - par.Ldom/2;
N = par.N;

tau_e = par.tau_e;  tau_i = par.tau_i;
a_ee = par.a_ee;  a_ei = par.a_ei;  a_ie = par.a_ie;  a_ii = par.a_ii;
th_e = par.th_e;  th_i = par.th_i;
sig_e = par.sig_e;
sig_i = par.sig_ratio*par.sig_e; 

F = @(u) 1./(1+exp(-u)); % sigmoid F
Kern = @(x,s) exp(-x.^2/s^2)/(s*sqrt(pi));  % kernel

Ke_hat = fft(ifftshift(Kern(x,sig_e)));
Ki_hat = fft(ifftshift(Kern(x,sig_i)));

% the jawn to be rootfound 
non_homog_syst = @(U) rhs_eqs(U, N, a_ee, a_ie, a_ei, a_ii, th_e, th_i, Ke_hat, Ki_hat, F, dx);

% both of these (spatially inhomog initial conditions) still lead to the
% spatially homog equilibrium
Ue0 = rand(1,par.N);  
Ui0 = rand(1,par.N); 
%Ue0 = 0.1*ones(1,N) + 1*exp(-x.^2/4); % bump
%Ui0 = 0.1*ones(1,N) + 1*exp(-x.^2/4);
U0 = [Ue0; Ui0]; 

[eqpt, fval] = fsolve(non_homog_syst, U0);

% define eq pt vals for each x for Ue, Ui
Ue_eq = eqpt(1:N);
Ui_eq = eqpt(N+1:end);

figure
plot(x, Ue_eq, 'b', 'LineWidth', 2);
hold on
plot(x, Ui_eq, 'r', 'LineWidth', 2)
xlabel('x')
ylabel('Equilibrium Value U(x)')
title('Spatial Equilibrium of Unforced Wilson-Cowan System')
legend('Excitatory (U_e)', 'Inhibitory (U_i)')


% finding eq for spatially homogeneous case: uncomment stuff below
% Ke_homog = 1/2*(erf(x(end)/sig_e) + erf(-x(1)/sig_e)); % spatially  homog, eq - kernels evaluate to constants (basically 1)
% Ki_homog = 1/2*(erf(x(end)/sig_i) + erf(-x(1)/sig_i));
% % 
% rhs_eq_syst = @(U) [
%     -U(1)/tau_e + F(a_ee*Ke_homog*U(1) - a_ie*Ki_homog*U(2)-th_e)/tau_e; 
%     -U(2)/tau_i + F(a_ei*Ke_homog*U(1) - a_ii*Ki_homog*U(2)-th_i)/tau_i;
%     ];
% 
%  U0 = [0.5; 0.25]; %  basically equivalent to enforcing spatially uniform
%  initial condition

% [eqpt, f] = fsolve(rhs_eq_syst, U0);

% we have the eq point!
% Ue = eqpt(1);
% Ui = eqpt(2);

% fprintf('eq point: (%.4f, %.4f), found near (%.4f, %.4f)\n', Ue, Ui, U0(1), U0(2))

% stability testing for the spatially homogeneous system
% a = (a_ee*Ue*(1-Ue) - 1)/tau_e;
% b = (-a_ie*Ue*(1-Ue))/tau_e;
% c = (a_ei*Ui*(1-Ui))/tau_i;
% d = (-a_ii*Ui*(1-Ui)-1)/tau_i;
% 
% p = [1, a+c, a*d-b*c]; % Represents x^2 - 2x + 1 = 0
% r = roots(p); 


% func for finding the eq without assuming spatial homogeneity:
function rhs_eqs_nonhomog = rhs_eqs(U, N, a_ee, a_ie, a_ei, a_ii, th_e, th_i, Ke_hat, Ki_hat, F, dx)
    Ue = U(1:N);
    Ui = U(N+1:end);
    size(Ki_hat)
    size(Ui)
    size(Ue)
    
    % use same code from model call for dealing with the conv.
    KeUe = dx*real(ifft(fft(Ue).*Ke_hat));
    KiUi = dx*real(ifft(fft(Ui).*Ki_hat));

    % define the odes
    Ue_de = -Ue + F(a_ee*KeUe - a_ie*KiUi - th_e);
    Ui_de = -Ui + F(a_ei*KeUe - a_ii*KiUi - th_i);
    
    rhs_eqs_nonhomog = [Ue_de; Ui_de];
end