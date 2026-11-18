function par = default_par()
    % spatial
    par.Ldom = 50; % length of domain
    par.N = 128; % number of points

    % time (ms)
    par.dt = 0.1; % timestep spacing
    par.tEnd = 2000; % final time

    % time constants
    par.tau_e = 10;
    par.tau_i = 20;

    % connection strengths
    par.a_ee = 10;
    par.a_ei = 12;
    par.a_ie = 8.5;
    par.a_ii = 3;

    % thresholds
    par.th_e = 2; % theta_e
    par.th_i = 3.5; % theta_i

    % kernel widths
    par.sig_e = 2;
    par.sig_ratio = 2.5; % sig_i = sig_ratio*sig_e

    % misc
    par.verbose = false;   % print kernel sanity checks
end