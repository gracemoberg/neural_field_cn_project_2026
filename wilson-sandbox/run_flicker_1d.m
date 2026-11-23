function [x,t,hist_e,hist_i] = run_flicker_1d(forcing_e, forcing_i, par, Ue, Ui)
    % access params from par struct
    N = par.N;  dt = par.dt;
    tau_e = par.tau_e;  tau_i = par.tau_i;
    a_ee = par.a_ee;  a_ei = par.a_ei;  a_ie = par.a_ie;  a_ii = par.a_ii;
    th_e = par.th_e;  th_i = par.th_i;
    sig_e = par.sig_e;
    sig_i = par.sig_ratio*par.sig_e;

    % spatial initialization
    dx = par.Ldom/N; % spacing
    x = (0:N-1)*dx - par.Ldom/2; % x vector, shift over so it is centered at zero

    % time (ms)
    t = 0:dt:par.tEnd; % time vector
    M = numel(t); % total number of time points

    % check the initial conditions match the grid
    assert(isequal(size(Ue), [1 N]) && isequal(size(Ui), [1 N]), 'Ue0 and Ui0 must be 1 x N row vectors (N = %d)', N);

    % inline functions
    F = @(u) 1./(1+exp(-u)); % sigmoid F
    Kern = @(x,s) exp(-x.^2/s^2)/(s*sqrt(pi));  % kernel

    % kernels - same each time so only need to compute once
    % use ifftshift to move the peak of the kernel from index N/2+1 to index 1
    % compute fft to transform kernel into Fourier space
    Ke_hat = fft(ifftshift(Kern(x,sig_e)));
    Ki_hat = fft(ifftshift(Kern(x,sig_i)));

    % verify kernel normalization and zero-frequency response
    if par.verbose
        fprintf('kernel norm (≈1): %.4f\n', sum(Kern(x,sig_e))*dx);
        fprintf('uniform gain (≈1): %.4f\n', real(Ke_hat(1))*dx);
        fprintf('max imag (≈0): %.2e\n', max(abs(imag(Ke_hat))));
    end

    % store Ue and Ui
    hist_e = zeros(M,N);
    hist_i = zeros(M,N);

    % use fft to transform Ue and Ui
    % shifted kernel real and even so fft is real
    % apply real() to remove roundoff errors
    % compute Ue and Ui using FW Euler
    for k = 1:M
        Ie = forcing_e(t(k), x);   % forcing at current time, amplitude already included
        Ii = forcing_i(t(k), x);
        KeUe = dx*real(ifft(fft(Ue).*Ke_hat));
        KiUi = dx*real(ifft(fft(Ui).*Ki_hat));
        Ue = Ue + dt/tau_e*(-Ue + F(a_ee*KeUe - a_ie*KiUi - th_e + Ie));
        Ui = Ui + dt/tau_i*(-Ui + F(a_ei*KeUe - a_ii*KiUi - th_i + Ii));
        hist_e(k,:) = Ue;
        hist_i(k,:) = Ui;
    end
end