%{
    ----------------------------------------------------------------------
    Author(s):    [Yannick DE BRUIJN, Michael FLOATER, Erik HILTUNEN]
    Date:         [August 2026]
    Description:  [First order phase correction for the GBZ preconditioner]
    ----------------------------------------------------------------------
%} 

clc
clear;
close all;


% ==== Parameters ====

    % --- Symbol ---
    m = 7;              % Truncation size for a_k
    p = 3.5;            % Decay rate upwards
    q = 4.8;            % Decay rate downwards

    % --- Discretisation ---
    num_lambda = 60;   % Number of sampling points on the band
    F_range    = 30;    % Fourier modes kept for the symmetrised symbol
    NF2        = 2^12;  % Grid for the correction symbol mu
    n_eps      = 120;   % Nodes (uniform in theta) where eps is evaluated
    collar     = 0.04;  % Edge collar excluded from the phase fit

    % --- Determinant recurrence fit ---
    n_phase = 60;       % Base size for the determinant recurrence fit
    n_seq   = 18;       % Number of consecutive sizes in the fit

    % --- Experiments ---
    DimList = [12 16 24 32 48 64 96];   % Sweep for d_sigma
    n_sign  = 48;       % Size used to fix the sign of the correction
    n_emp   = 96;       % Size for the empirical cross-check
    fs      = 18;       % Fontsize for annotation


% ==== Generate m-banded Dummy Toeplitz Matrix ====

    col = zeros(m,1);
    row = zeros(1,m+1);

    % --- Add noise to the coeffiients ---
    ai = 1;
    bi = 1;

    col(1) = 1;
    row(1) = 1;

    % --- Populate above and below diagonals ---
    for k = 1:m
        r = ai + (bi-ai)*rand;
        col(k) = r / ((k+1)^q);
        r = ai + (bi-ai)*rand;
        row(k+1) = r / ((k+1)^p);
    end

    % --- Coefficients of the symbol function ---
    a = [ col(end:-1:1)', row];


% ==== Approximate the open limit ====

    [lambda_interval(1), lambda_interval(2)] = open_limit(a);


% ==== Approximate conjugate root set Lambda(f) ====

    % --- Chebyshev clustering: the roots collide like sqrt at the edges ---
    t = linspace(0, 1, num_lambda);
    lambda_vals = lambda_interval(1) + (lambda_interval(2) - lambda_interval(1)) ...
                  * (0.5 * (1 - cos(pi*t)));


% ==== Compute discrete points on GBZ ====

    % --- Roots of z^m (f(z) - lambda), sorted by ascending modulus ---
    P_coeffs_matrix = repmat(a(:), 1, num_lambda);
    P_coeffs_matrix(m + 1, :) = P_coeffs_matrix(m + 1, :) - lambda_vals;

    all_roots = zeros(2*m, num_lambda);
    for k = 1:num_lambda
        all_roots(:, k) = roots(P_coeffs_matrix(:, k));
    end

    [~, sort_idx] = sort(abs(all_roots), 1);
    linear_idx = sort_idx + (0:num_lambda-1) * (2*m);
    all_roots_sorted = all_roots(linear_idx);

    % --- Keep the central pair, and only where the moduli agree ---
    candidate_roots = all_roots_sorted(m:m+1, :);

    tolerance = 1e-8;  % keep in mind rootsolver has its limits
    mod_n  = abs(candidate_roots(1, :));
    mod_n1 = abs(candidate_roots(2, :));
    similar_modulus = abs(mod_n - mod_n1) ./ max(mod_n, mod_n1) < tolerance;


% ==== Band function lambda(theta) and its inverse ====

    % --- Table (theta, lambda) read off from the central pair ---
    keep  = find(similar_modulus);
    th_g  = zeros(1, numel(keep));
    lam_g = lambda_vals(keep);
    for i = 1:numel(keep)
        zpair = candidate_roots(:, keep(i));
        [~, iP] = max(imag(zpair));
        th_g(i) = abs(angle(zpair(iP)));
    end

    [th_g, ord] = sort(th_g);
    lam_g = lam_g(ord);
    kp = [true, diff(th_g) > 1e-10];
    th_g = th_g(kp);  lam_g = lam_g(kp);

    % --- Pin the edges, where theta = 0 and theta = pi ---
    th_b  = [0, th_g, pi];
    lam_b = [lambda_interval(1), lam_g, lambda_interval(2)];
    if lam_g(1) > lam_g(end)
        lam_b = [lambda_interval(2), lam_g, lambda_interval(1)];
    end
    lamS = griddedInterpolant(th_b, lam_b, 'spline');

    % --- Monotone inverse theta(lambda), built once and reused below.
    %     unique() drops the duplicated nodes and sorts ascending, so the
    %     pairing with th_b(iu) stays valid in either orientation ---
    [lam_inv, iu] = unique(lam_b(:));
    th_inv = th_b(iu);  th_inv = th_inv(:);
    thOf = @(x) interp1(lam_inv, th_inv, min(max(x, lam_inv(1)), lam_inv(end)), 'linear', 'extrap');


% ==== Resample the GBZ ====

    candidate_roots = candidate_roots(:, similar_modulus);
    openLimit = [candidate_roots(1, :), candidate_roots(2, :)];

    phase = angle(openLimit(:));
    [~, sortIdx] = sort(phase);
    openLimit_sorted = openLimit(sortIdx);

    openLimit_sorted = merge_close_points(openLimit_sorted, 1e-1);

    nols  = length(openLimit_sorted);
    nolsh = nols/2;
    olsh  = openLimit_sorted(nolsh+2:nols);
    mnew  = num_lambda*5;
    mnew2 = mnew/2+1;

    interp_openLimit = resample2(olsh, mnew2);

    % --- Shift resample2 output to standard FFT order 0..2pi ---
    interp_openLimit = circshift(interp_openLimit, -(mnew2-2));


% ==== Compute Fourier coefficients of f(p(z)) numerically ====

    k_values = -m:m;
    powers_matrix = interp_openLimit(:).^(-k_values);
    fp_values = powers_matrix * a(:);
    fp_values = real(fp_values);

    FourierFP = fourier_coefficients_spectral(fp_values, F_range);
    gco = FourierFP(F_range+1:end);          % one sided, gco(k+1) = g_k


% ==== Phase profile from determinant data only ====
    % For each lambda we fit the two term recurrence satisfied by the
    % determinant sequence and read off the boundary phase. Both sides
    % share the same theta(lambda) (same band function), so the whole
    % first order discrepancy is phi_g - phi_f.

    [Sidx, Cidx] = widom_subsets(m);

    th_raw = linspace(collar, pi - collar, n_eps);
    lam_e  = lamS(th_raw);
    eps_r  = zeros(size(th_raw));
    th_chk = zeros(size(th_raw));

    nvals = n_phase + (0:n_seq-1);

    for i = 1:numel(th_raw)
        lam = lam_e(i);

        % --- f side: scaled Widom determinant, D_n / G^n ---
        vf = zeros(1, n_seq);
        for j = 1:n_seq
            vf(j) = widom_det(a, lam, nvals(j), Sidx, Cidx);
        end
        [th_f, phi_f] = phase_from_seq(vf, nvals);

        % --- g side: LU log determinant of the symmetric banded matrix ---
        vg = det_sequence_sym(gco, lam, nvals);
        [th_h, phi_g] = phase_from_seq(vg, nvals);

        % --- The two thetas must agree: this is Theorem 2.10 numerically ---
        th_chk(i) = abs(th_f - th_h);

        % --- Phase difference, principal value (zeros are mod pi) ---
        d = phi_g - phi_f;
        d = mod(d + pi/2, pi) - pi/2;
        eps_r(i) = d;
    end

    fprintf('max |theta_f - theta_g| over the band: %.2e  (should be ~0)\n', ...
            max(th_chk));
    fprintf('max |eps| = %.3e rad\n', max(abs(eps_r)));

    % --- Unwrap along the band and pin the edges, where eps vanishes ---
    eps_r = unwrap_small(eps_r);
    th_e  = [0, th_raw, pi];
    eps_e = [0, eps_r, 0];


% ==== First order correction term ====

    th_full  = 2*pi*(0:NF2-1)'/NF2;
    th_fold  = pi - abs(pi - th_full);
    sgn_full = sign(sin(th_full) + (sin(th_full) == 0));

    % --- lambda'(theta) by spectral differentiation of the g_k ---
    kk = 1:F_range;
    dlam_full = -2*sin(th_full*kk) * (kk(:) .* gco(2:F_range+1));

    % --- eps on the folded grid (smooth now, pchip is safe) ---
    eps_fold = interp1(th_e, eps_e, th_fold, 'pchip');
    eps_full = sgn_full .* eps_fold;

    % --- Correction symbol mu = lambda' * eps (odd * odd = even) ---
    mu_full    = dlam_full .* eps_full;
    Fourier_mu = fourier_coefficients_spectral(mu_full, 2*F_range);


% ==== Fix the sign with one test at moderate n ====

    eig_f = widom_spectrum(a, n_sign, lambda_interval, Sidx, Cidx);
    T_b   = fourier_to_toeplitz(FourierFP, n_sign);
    T_mu  = fourier_to_toeplitz(Fourier_mu, n_sign);

    d_plus  = sum(abs(eig_f(:) - sort(eig(T_b + T_mu/n_sign))));
    d_minus = sum(abs(eig_f(:) - sort(eig(T_b - T_mu/n_sign))));

    if d_minus < d_plus
        Fourier_mu = -Fourier_mu;
        eps_e      = -eps_e;
    end
    fprintf('sign test at n=%d:  d(+) = %.3e,  d(-) = %.3e\n', ...
            n_sign, d_plus, d_minus);


% ==== Predicted plateau of the uncorrected distance ====

    lam_e2 = lamS(th_e);
    c_inf  = trapz(lam_e2, abs(eps_e)) / pi;
    fprintf('predicted plateau  c_inf = %.5e\n\n', abs(c_inf));


% ==== d_sigma sweep ====

    d_unc = zeros(size(DimList));
    d_cor = zeros(size(DimList));

    fprintf('    n     d_sigma(T,Tb)    d_sigma(T,Tb + mu/n)   ratio\n');
    for iq = 1:length(DimList)
        DimT = DimList(iq);

        % --- Reference spectrum of the non-Hermitian matrix ---
        eig_f = widom_spectrum(a, DimT, lambda_interval, Sidx, Cidx);

        % --- Symmetrised matrix, without and with the correction ---
        T_b  = fourier_to_toeplitz(FourierFP, DimT);
        T_mu = fourier_to_toeplitz(Fourier_mu, DimT);

        eig_b = sort(eig(T_b));
        eig_c = sort(eig(T_b + T_mu/DimT));

        d_unc(iq) = sum(abs(eig_f(:) - eig_b(:)));
        d_cor(iq) = sum(abs(eig_f(:) - eig_c(:)));

        fprintf('  %4d   %13.5e    %13.5e    %8.2f\n', ...
                DimT, d_unc(iq), d_cor(iq), d_unc(iq)/d_cor(iq));
    end


% ==== Empirical cross-check of eps ====
    % The measured displacement of the quantisation grid, rescaled by n,
    % must reproduce the phase difference computed above.

    eig_f = widom_spectrum(a, n_emp, lambda_interval, Sidx, Cidx);
    eig_b = sort(eig(fourier_to_toeplitz(FourierFP, n_emp)));

    sf = thOf(eig_f(:));
    sH = thOf(eig_b(:));
    eps_emp = (n_emp + 1) * (sf - sH);

    %{
    figure;
    hold on;
    plot(sH, eps_emp, 'r.', 'MarkerSize', 9);
    plot(th_e, eps_e, 'b-', 'LineWidth', 2);
    grid on;
    box on;
    xlabel('$\theta$', 'Interpreter', 'latex', 'FontSize', fs);
    ylabel('$\varepsilon(\theta)$', 'Interpreter', 'latex', 'FontSize', fs);
    legend({'$(n{+}1)(s^f_j - s^H_j)$', '$\varphi_g - \varphi_f$'}, ...
           'Interpreter', 'latex', 'FontSize', fs-4, 'Location', 'best');
    set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', fs);
    set(gcf, 'Position', [100, 100, 500, 300]);
    hold off;

    %}


% ==== Plot the l1 distance before and after the correction ====

    figure;
    loglog(DimList, d_unc, 'k-o', 'LineWidth', 1.5);
    hold on;
    %loglog(DimList, d_cor, 'b-s', 'LineWidth', 1.5);
    loglog(DimList, abs(c_inf)*ones(size(DimList)), 'g:', 'LineWidth', 1.5);
    %loglog(DimList, d_cor(1)*DimList(1)./DimList, 'r--', 'LineWidth', 1.5);
    grid on;
    box on;
    xlabel('Matrix size $n$', 'Interpreter', 'latex', 'FontSize', fs);
    ylabel('$d_\sigma$', 'Interpreter', 'latex', 'FontSize', fs);
    legend({'$d_\sigma(\mathbf{T}_n(f), \mathbf{T}_n(f\circ p))$', ...
            %'$d_\sigma(\mathbf{T}_n(f), \mathbf{T}_n(f\circ p) + \frac{1}{n}\mathbf{T}_n(\lambda''\varepsilon))$', ...
            '$c_\infty$ prediction', ...
            %'$\mathcal{O}(n^{-1})$'
            }, ...
            'Interpreter', 'latex', 'Location', 'southwest');
    set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', fs);
    set(gcf, 'Position', [100, 100, 500, 300]);
    hold off;


%% --- Defining functions ---


function [th, phi] = phase_from_seq(v, nvals)
% Boundary phase from a determinant sequence. On the band
%     v_n = A G^n cos(n th + phi),
% so v satisfies v_{n+1} = p v_n + q v_{n-1} with characteristic roots
% G exp(+-i th). We fit p,q by least squares, then the complex amplitude.
% Any effective size offset (n + kappa) is absorbed into phi, which is
% exactly what the correction needs.

    v = v(:);
    K = numel(v);

    % --- Fit the two term recurrence ---
    M   = [v(2:K-1), v(1:K-2)];
    rhs = v(3:K);
    pq  = M \ rhs;

    mu = roots([1, -pq(1), -pq(2)]);
    if abs(imag(mu(1))) < 1e-12
        th = NaN;  phi = NaN;  return;    % no oscillation: off the band
    end
    G  = abs(mu(1));
    th = abs(angle(mu(1)));

    % --- Fit the complex amplitude: v_k = Re(C mu^k) ---
    kk = (0:K-1)';
    Bs = [G.^kk .* cos(kk*th), -G.^kk .* sin(kk*th)];
    C  = Bs \ v;

    % --- Phase at index n = nvals(1) + k, so unwind the base offset ---
    psi0 = atan2(C(2), C(1));
    phi  = psi0 - nvals(1)*th;
    phi  = mod(phi + pi, 2*pi) - pi;
end


function v = det_sequence_sym(gco, lam, nvals)
% Determinants det(T_n(g) - lam) for consecutive n, scaled by the first
% one so the sequence is O(1). The matrix is symmetric and banded, hence
% well conditioned, so a plain LU log determinant is enough: no roots, no
% truncated colleague matrix, no root selection.

    K  = numel(nvals);
    ld = zeros(1, K);
    sg = zeros(1, K);

    for j = 1:K
        n = nvals(j);

        % --- Symmetric banded Toeplitz matrix at this size ---
        c = zeros(n, 1);
        L = min(numel(gco), n);
        c(1:L) = gco(1:L);
        A = toeplitz(c, c) - lam*eye(n);

        % --- Log determinant, sign kept separately ---
        [~, U, pvec] = lu(A, 'vector');
        du = diag(U);
        ld(j) = sum(log(abs(du)));
        sg(j) = prod(sign(du)) * perm_sign(pvec);
    end

    v = sg .* exp(ld - ld(1));
end


function s = perm_sign(pvec)
% Sign of a permutation given in vector form.

    n = numel(pvec);
    seen = false(1, n);
    s = 1;
    for i = 1:n
        if ~seen(i)
            j = i;
            len = 0;
            while ~seen(j)
                seen(j) = true;
                j = pvec(j);
                len = len + 1;
            end
            if mod(len - 1, 2) == 1
                s = -s;
            end
        end
    end
end


function y = unwrap_small(x)
% Remove the occasional pi jump in a phase difference that is known to be
% small and continuous along the band.

    y = x;
    for i = 2:numel(y)
        while y(i) - y(i-1) >  pi/2, y(i) = y(i) - pi; end
        while y(i) - y(i-1) < -pi/2, y(i) = y(i) + pi; end
    end
end


function ck = fourier_coefficients_spectral(fp_values, K)
% Two sided Fourier coefficients of an even real sample vector given in
% standard FFT order, truncated at K modes.

    M = length(fp_values);

    % --- Use standard (uniform) FFT ---
    fft_result = fft(fp_values) / M;

    % --- Extract coefficients ---
    ck = zeros(2*K+1, 1);
    ck(K+1)     = fft_result(1);
    ck(K+2:end) = fft_result(2:K+1);
    ck(1:K)     = fft_result(end-K+1:end);

    % --- Even real samples in standard FFT order => real coefficients ---
    ck = real(ck);
end


function T = fourier_to_toeplitz(a, dimT)
% Finite section T_n(h) from the two sided coefficients of h.

    K = (length(a) - 1) / 2;
    a_0 = a(K+1);

    col = zeros(dimT, 1);
    col(1) = a_0;
    row = zeros(1, dimT);
    row(1) = a_0;

    for k = 1:min(K, dimT-1)
        col(k+1) = a(K+1-k);
        row(k+1) = a(K+1+k);
    end

    T = toeplitz(col, row);
end


function [lambdaL, lambdaR] = open_limit(a)
% Endpoints of the open limit by bisection on the modulus gap
% |z_{m+1}(lambda)| - |z_m(lambda)|, bracketed by f(+-1) and the spectrum
% of a small finite section.

    tol      = 1e-8;   % criterion: how flat is "zero"
    min_step = 1e-5;   % geometric: stop bisecting here (<< tol/|F'|)

    m = (length(a)-1)/2;

    % --- Coefficients of Q in ascending powers ---
    coeff_Q = zeros(1, 2*m+1);
    for j = -m:m
        coeff_Q(m-j+1) = a(j+m+1);
    end

    % --- Outer brackets: the edge of the winding region ---
    f1 = sum(a);
    j = -m:m;
    f_minus1 = sum(a .* (-1).^j);

    % --- Inner brackets: a small finite section ---
    T = fourier_to_toeplitz(a, 10);
    eigT = sort(eig(T));

    fmin = min(eigT);
    fmax = max(eigT);

    function F = diff_mod(lambda)
        c = coeff_Q;
        c(m+1) = c(m+1) - lambda;
        r = roots(c);
        r = sort(abs(r));
        F = abs(r(m+1)) - abs(r(m)) - abs(imag(r(m+1)));
    end

    function inside = in_interval(lambda)
        inside = abs(diff_mod(lambda)) <= tol;
    end

    function endpoint = find_endpoint(inside, outside)
        while abs(outside - inside) > min_step
            mid = (inside + outside) / 2;
            if in_interval(mid)
                inside = mid;
            else
                outside = mid;
            end
        end
        endpoint = inside;
    end

    lambdaR = find_endpoint(fmax, f1);
    lambdaL = find_endpoint(fmin, f_minus1);
end


function merged = merge_close_points(openLimit, tol)
% Merge the duplicated endpoints produced when the two branches of the
% central pair are concatenated.

    n = length(openLimit);
    if n == 0
        merged = [];
        return;
    end

    merged = openLimit;

    % --- Wrap around duplicate ---
    if abs(openLimit(end) - openLimit(1)) <= tol
        merged(1) = mean([openLimit(1), openLimit(end)]);
        merged(end) = [];
    end

    % --- Duplicate where the two branches meet ---
    if n >= 2
        mid1 = floor(n / 2);
        mid2 = mid1 + 1;

        if abs(merged(mid2) - merged(mid1)) <= tol
            merged(mid1) = mean([merged(mid1), merged(mid2)]);
            merged(mid2) = [];
        end
    end
end


function [Sidx, Cidx] = widom_subsets(m)
% All splittings of the 2m roots into m "inside" and m "outside".

    Sidx = nchoosek(1:2*m, m);
    Cidx = zeros(size(Sidx));
    for i = 1:size(Sidx,1)
        Cidx(i,:) = setdiff(1:2*m, Sidx(i,:));
    end
end


function v = widom_det(a, lambda, n, Sidx, Cidx)
% Scaled Widom determinant of T_n(f) - lambda, i.e. D_n / G^n. All terms
% are computed in log scale relative to the dominant modulus, so the
% result is O(1) and real on the band.

    m = (length(a)-1)/2;

    % --- Roots of z^m (f(z) - lambda), ascending modulus ---
    c = a;
    c(m+1) = c(m+1) - lambda;
    z = roots(c);
    [~, o] = sort(abs(z));
    z = z(o);

    lz   = log(z);
    lead = log((-1)^m * a(1));          % a(1) = a_{-m}

    % --- One term per splitting, weight and prefactor in log scale ---
    zS  = z(Sidx);
    zC  = z(Cidx);
    slS = sum(lz(Sidx), 2);

    logw = lead + slS;
    D    = reshape(zS, [], m, 1) - reshape(zC, [], 1, m);
    lc   = m*slS - sum(sum(log(D), 3), 2);

    v = real(sum(exp(lc + n*(logw - max(real(logw))))));
end


function ev = widom_spectrum(a, n, lambda_interval, Sidx, Cidx)
% Exact spectrum of T_n(f) via sign changes of the Widom determinant.

    oversamp = 24;
    pad  = 0.02*(lambda_interval(2) - lambda_interval(1));
    grid = linspace(lambda_interval(1) - pad, lambda_interval(2) + pad, oversamp*n);

    % --- Bracket every sign change on the oversampled grid ---
    v = arrayfun(@(L) widom_det(a, L, n, Sidx, Cidx), grid);
    s = sign(v);
    ix = find(s(1:end-1).*s(2:end) < 0);

    lo = grid(ix);
    hi = grid(ix+1);
    flo = v(ix);

    % --- Vectorised bisection on all brackets at once ---
    for it = 1:48
        mid = 0.5*(lo + hi);
        fm  = arrayfun(@(L) widom_det(a, L, n, Sidx, Cidx), mid);
        left = flo.*fm < 0;
        hi(left)   = mid(left);
        lo(~left)  = mid(~left);
        flo(~left) = fm(~left);
    end

    ev = sort(0.5*(lo + hi))';
end