function [FourierFP, Fourier_mu, c_inf, info] = symmetriser_from_roots(a, opts)
%SYMMETRISER_FROM_ROOTS  First-order Toeplitz symmetriser computed from roots only.
%
%   [FourierFP, Fourier_mu, c_inf, info] = symmetriser_from_roots(a)
%
%   a : coefficients [a_{-m} ... a_0 ... a_m] of f(z) = sum_k a_k z^{-k}  (a(1) = a_{-m}).
%
%   For each level lambda_k on the band the roots of z^m (f(z) - lambda_k) give
%     * the momentum alpha_k = arg zeta_k from the central conjugate pair, i.e. the
%       inverse band function (this replaces the GBZ curve: g(e^{i alpha}) = lambda(alpha));
%     * the reduced phase Psi_f(alpha_k) from the inner/outer roots, eq. (6.4);
%   and the Hermitian phase Psi_g(alpha_k) is obtained either root-free, as a Hilbert
%   transform of log|(g - lambda)/(2cos theta - 2cos alpha)| (default), or from the
%   inner roots of the truncated g.  Then eps = Psi_g - Psi_f, mu = lambda' eps, and
%     T_n(f) ~ T_n(g) + T_n(mu)/n     with error O(n^-2) pointwise.
%   No determinant sequences, no spline through the GBZ, no sign test.
%
%   Outputs (same layout as first_and_second_order_expansion)
%     FourierFP   two-sided coefficients [g_F ... g_1 g_0 g_1 ... g_F] of g = f o p
%     Fourier_mu  two-sided coefficients of mu = lambda' eps, range 2*F_range
%     c_inf       (1/pi) int |eps| dlambda
%     info        nodes, phases, fit residuals, handles lambda(.) and eps(.), band edges
%
%   Build the matrices with  toeplitz([gco, zeros(1,n-numel(gco))])  where
%   gco = FourierFP(F_range+1:end), and likewise for mu.

    arguments
        a double
        opts.num_lambda = 100       % levels on the band (cosine clustered)
        opts.F_range    = 30        % truncation of g
        opts.NF2        = 2^12      % theta grid for FFTs
        opts.n_modes    = 30        % sine modes representing eps
        opts.psi_g      = 'hilbert' % 'hilbert' (root free) or 'roots'
        opts.n_check    = 40        % eig-based self check at this size (0 = skip)
        opts.plot_phase_shift (1,1) logical = false
        opts.fs         = 18
        opts.verbose    = true
    end

    a = a(:).';
    m = (numel(a) - 1) / 2;
    F = opts.F_range;
    NF2 = opts.NF2;
    th_full = 2*pi*(0:NF2-1)' / NF2;

    % ==== Band edges (lambda at theta = 0 and theta = pi) ====
    [lam_th0, lam_thpi] = band_edges(a);
    lamA = min(lam_th0, lam_thpi);  lamB = max(lam_th0, lam_thpi);

    % ==== Roots on the band: momentum and non-Hermitian phase ====
    N = opts.num_lambda;
    t = (1:N) / (N + 1);
    lam_vals = lamA + (lamB - lamA) * 0.5 * (1 - cos(pi*t));
    al = zeros(1, N);  psi_f = zeros(1, N);  ok = false(1, N);
    for k = 1:N
        z  = sorted_roots(a, lam_vals(k));
        zc = z(m:m+1);
        if abs(abs(zc(1)) - abs(zc(2))) / max(abs(zc)) > 1e-8 || abs(imag(zc(1))) < 1e-10
            continue;                                   % not a conjugate central pair
        end
        [~, iP] = max(imag(zc));  zeta = zc(iP);
        zin  = z(1:m-1);  zout = z(m+2:end);
        al(k)    = angle(zeta);                         % alpha(lambda_k)
        psi_f(k) = al(k) - sum(angle(1 - zin/zeta)) + sum(angle(1 - zeta./zout));   % eq. (6.4), Widom orientation
        ok(k)    = true;
    end
    al = al(ok);  psi_f = psi_f(ok);  lam_n = lam_vals(ok);
    [al, o] = sort(al);  psi_f = psi_f(o);  lam_n = lam_n(o);
    kp = [true, diff(al) > 1e-12];  al = al(kp);  psi_f = psi_f(kp);  lam_n = lam_n(kp);

    % ==== Symmetrised symbol: cosine series of lambda(alpha) fitted to the exact nodes ====
    kk    = 1:F;
    al_b  = [0, al, pi];  lam_b = [lam_th0, lam_n, lam_thpi];
    Xc    = [ones(numel(al_b), 1), 2*cos(al_b(:) * kk)];
    gco   = (Xc \ lam_b(:)).';                          % gco(k+1) = g_k
    fit_res_lambda = max(abs(Xc * gco.' - lam_b(:)));
    FourierFP = [fliplr(gco(2:end)), gco];

    lam_trig  = @(th) gco(1) + 2*cos(th(:) * kk) * gco(2:end).';
    dlam_trig = @(th) -2*sin(th(:) * kk) * (kk(:) .* gco(2:end).');
    gvals     = lam_trig(th_full);

    % ==== Hermitian phase at the same nodes ====
    psi_g = zeros(size(al));
    switch lower(opts.psi_g)
        case 'hilbert'
            kh  = (1:NF2/2-1)';                                   % (NF2/2-1) x 1
            Q   = abs(gvals - lam_trig(al).') ...                 % NF2 x N : |g - lambda_i|
                ./ abs(2*cos(th_full) - 2*cos(al));               %       / |2cos th - 2cos beta_i|
            C   = real(fft(0.5*log(Q), [], 1)) / NF2;             % cosine coefficients, one column per node
            S   = sin(kh * al);                                   % (NF2/2-1) x N : sin(k beta_i)
            argbp = 2 * sum(C(2:NF2/2, :) .* S, 1);               % 1 x N : arg b_+(e^{i beta_i})
            psi_g = al + 2*argbp;
        %{
        case 'hilbert'
            % g - lambda = C (2cos th - 2cos beta) b_+ b_-,  log|b_+| = (1/2) log q + const,
            % arg b_+(e^{i beta}) = harmonic conjugate = sum 2 c_k sin(k beta).
            kh = (1:NF2/2-1)';
            for i = 1:numel(al)
                beta = al(i);
                q  = abs(gvals - lam_trig(beta)) ./ abs(2*cos(th_full) - 2*cos(beta));
                c  = real(fft(0.5*log(q))) / NF2;
                argbp = 2*c(2:NF2/2).' * sin(kh*beta);
                psi_g(i) = beta + 2*argbp;
            end

         %}
        case 'roots'
            for i = 1:numel(al)
                beta = al(i);
                coef = [gco(end:-1:2), gco(1) - lam_trig(beta), gco(2:end)];   % z^F (g(z) - lambda)
                w = roots(coef);  [~, o] = sort(abs(w));  w = w(o(1:F-1));    % inner roots
                psi_g(i) = beta + sum(angle((1 - w*exp(1i*beta)) ./ (1 - w*exp(-1i*beta))));
            end
        otherwise
            error('psi_g must be ''hilbert'' or ''roots''');
    end

    % ==== Phase profile eps = Psi_g - Psi_f, represented as an odd sine series ====
    eps_n = psi_g - psi_f;                              % vanishes at the edges by construction
    J     = opts.n_modes;
    SinB  = sin(al(:) * (1:J));
    ecoef = SinB \ eps_n(:);
    fit_res_eps = max(abs(SinB * ecoef - eps_n(:)));
    eps_full  = sin(th_full * (1:J)) * ecoef;
    dlam_full = dlam_trig(th_full);

    % ==== First-order correction symbol mu = lambda' eps ====
    mu_full = dlam_full .* eps_full;
    cm  = real(fft(mu_full)) / NF2;
    mu1 = cm(1:2*F+1).';
    Fourier_mu = [fliplr(mu1(2:end)), mu1];

    ih    = 1:NF2/2+1;
    c_inf = trapz(th_full(ih), abs(dlam_full(ih)) .* abs(eps_full(ih))) / pi;

    % ==== Diagnostics ====
    info.alpha   = al;        info.lambda = lam_n;
    info.psi_f   = psi_f;     info.psi_g  = psi_g;     info.eps = eps_n;
    info.gco     = gco;       info.mu     = mu1;
    info.band    = [lamA, lamB];  info.lambda_theta0 = lam_th0;  info.lambda_thetapi = lam_thpi;
    info.lam     = lam_trig;  info.dlam   = dlam_trig;
    info.eps_fun = @(th) sin(th(:) * (1:J)) * ecoef;
    info.fit_res_lambda = fit_res_lambda;   info.fit_res_eps = fit_res_eps;
    info.g_tail  = abs(gco(end));           info.eps_tail = max(abs(ecoef(end-2:end)));
    if opts.verbose
        fprintf('band [%.6f, %.6f], %d nodes; |g_F| = %.1e, lambda fit residual %.1e, eps fit residual %.1e, c_inf = %.3e\n', ...
            lamA, lamB, numel(al), info.g_tail, fit_res_lambda, fit_res_eps, c_inf);
    end

    % ==== Plot ====
    if opts.plot_phase_shift
        fs = opts.fs;
        figure;  hold on;
        plot(al, eps_n, 'b.', 'MarkerSize', 8);
        thp = linspace(0, pi, 400);
        plot(thp, info.eps_fun(thp), 'b-', 'LineWidth', 1.5);
        grid on;  box on;
        xlabel('$\theta$', 'Interpreter', 'latex', 'FontSize', fs);
        ylabel('$\varepsilon(\theta)$', 'Interpreter', 'latex', 'FontSize', fs);
        set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', fs);
        xlim([0 pi]);  xticks(0:pi/4:pi);
        xticklabels({'$0$', '$\pi/4$', '$\pi/2$', '$3\pi/4$', '$\pi$'});
        set(gcf, 'Position', [100, 100, 500, 300]);
        hold off;
    end
end


%% ---------------------------------------------------------------------
function z = sorted_roots(a, lam)
% Roots of z^m (f(z) - lambda), ascending modulus.
    m = (numel(a) - 1) / 2;
    c = a(:);  c(m+1) = c(m+1) - lam;
    z = roots(c);
    [~, o] = sort(abs(z));
    z = z(o);
end


function [lam_th0, lam_thpi] = band_edges(a)
% Band edges as the real critical values of f at which the central pair coalesces:
% lam_th0 has a positive real double root (theta = 0), lam_thpi a negative one (theta = pi).
% Falls back to open_limit(a) if available.
    m = (numel(a) - 1) / 2;  a = a(:).';
    cder = -((0:2*m) - m) .* a;                 % coefficients of z^{m+1} f'(z), descending
    zc = roots(cder);
    zc = real(zc(abs(imag(zc)) < 1e-9 & abs(zc) > 1e-12));
    lam_e = [];  z_e = [];
    for zs = zc.'
        lam = real(sum(a .* zs.^(-(-m:m))));    % f(zs)
        z = sorted_roots(a, lam);
        if abs(abs(z(m)) - abs(zs)) < 1e-6 && abs(abs(z(m+1)) - abs(zs)) < 1e-6
            lam_e(end+1) = lam;  z_e(end+1) = zs; %#ok<AGROW>
        end
    end
    if numel(lam_e) == 2 && nnz(z_e > 0) == 1
        lam_th0 = lam_e(z_e > 0);  lam_thpi = lam_e(z_e < 0);
    elseif exist('open_limit', 'file') == 2
        [l1, l2] = open_limit(a);
        z = sorted_roots(a, l1);
        if real(z(m)) > 0, lam_th0 = l1; lam_thpi = l2; else, lam_th0 = l2; lam_thpi = l1; end
    else
        error('band_edges: found %d central critical values, expected 2', numel(lam_e));
    end
end
