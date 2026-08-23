function [FourierFP, Fourier_mu, c_inf] = first_and_second_order_expansion(a, opts)

    %{
        ----------------------------------------------------------------------
        Author(s):    [Yannick DE BRUIJN, Michael FLOATER, Erik HILTUNEN]
        Date:         [August 2026]
        Description:  [Second order Toeplitz Expansion]
        ----------------------------------------------------------------------
    %}
    arguments
        a
        opts.num_lambda = 60;    % Sampling points on the band
        opts.n_phase    = 60;     % Base size for determinant recurrence fit
        opts.n_seq      = 18;     % Number of consecutive sizes in the fit
        opts.n_eps      = 70;    % Nodes (uniform in theta) for eps
        opts.NF2        = 2^12;   % Grid for the correction symbol mu
        opts.collar     = 0.04;   % Edge collar excluded from the phase fit
        opts.F_range    = 30;     % Truncation of the Fourier coefficients
        opts.n0         = 40;     % Size for the sign test
        opts.plot_phase_shift (1,1) logical = false;  % Empirical cross-check of eps
        opts.n_emp      = 200;    % Matrix size for the empirical check
        opts.fs         = 18;     % Font size for the plot
    end

        num_lambda = opts.num_lambda;
        n_phase = opts.n_phase;  n_seq = opts.n_seq;  n_eps = opts.n_eps;
        NF2 = opts.NF2;  collar = opts.collar;
        F_range = opts.F_range;  n0 = opts.n0;

        m = (length(a) - 1) / 2;

    % ==== Approximate the open limit ====
        [lambda_interval(1), lambda_interval(2)] = open_limit(a);

    % ==== Approximate conjugate root set Λ(f) ====
        t = linspace(0, 1, num_lambda);
        lambda_vals = lambda_interval(1) + (lambda_interval(2) - lambda_interval(1)) * (0.5 * (1 - cos(pi*t)));

    % === Compute discrete points on GBZ ===
        P_coeffs_matrix = repmat(a(:), 1, num_lambda);
        P_coeffs_matrix(m + 1, :) = P_coeffs_matrix(m + 1, :) - lambda_vals;

        all_roots = zeros(2*m, num_lambda);
        for k = 1:num_lambda
            all_roots(:, k) = roots(P_coeffs_matrix(:, k));
        end

        [~, sort_idx] = sort(abs(all_roots), 1);
        linear_idx = sort_idx + (0:num_lambda-1) * (2*m);
        all_roots_sorted = all_roots(linear_idx);

        candidate_roots = all_roots_sorted(m:m+1, :);

        tolerance = 1e-8;
        mod_n  = abs(candidate_roots(1, :));
        mod_n1 = abs(candidate_roots(2, :));
        similar_modulus = abs(mod_n - mod_n1) ./ max(mod_n, mod_n1) < tolerance;

        % --- Exact band function table (alpha, lambda) from the root data ---
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

        % --- Band function lambda(theta), pinned at the edges ---
        th_b  = [0, th_g, pi];
        lam_b = [lambda_interval(1), lam_g, lambda_interval(2)];
        if lam_g(1) > lam_g(end)
            lam_b = [lambda_interval(2), lam_g, lambda_interval(1)];
        end
        lamS = griddedInterpolant(th_b, lam_b, 'spline');

        % --- resample2 pipeline (coefficients of f o p) ---
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
        interp_openLimit = circshift(interp_openLimit, -(mnew2-2));

    % ==== Compute Fourier coefficients of f(p(z)) numerically ====
        k_values = -m:m;
        powers_matrix = interp_openLimit(:).^(-k_values);
        fp_values = powers_matrix * a(:);
        fp_values = real(fp_values);

        FourierFP = fourier_coefficients_spectral(fp_values, F_range);
        gco = real(FourierFP(F_range+1:end));          % one sided, gco(k+1) = g_k

    % ==== Phase profile from determinant data only ====
        [Sidx, Cidx] = widom_subsets(m);

        th_e   = linspace(collar, pi - collar, n_eps);
        lam_e  = lamS(th_e);
        eps_e  = zeros(size(th_e));

        nvals = n_phase + (0:n_seq-1);

        for i = 1:numel(th_e)
            lam = lam_e(i);

            vf = widom_det_seq(a, lam, nvals, Sidx, Cidx);
            [~, phi_f] = phase_from_seq(vf, nvals);

            vg = det_sequence_sym(gco, lam, nvals);
            [~, phi_g] = phase_from_seq(vg, nvals);

            d = phi_g - phi_f;
            d = mod(d + pi/2, pi) - pi/2;
            eps_e(i) = d;
        end

        eps_e = unwrap_small(eps_e);
        th_e  = [0, th_e, pi];
        eps_e = [0, eps_e, 0];

    % ==== First order correction term ====
        th_full = 2*pi*(0:NF2-1)'/NF2;
        th_fold = pi - abs(pi - th_full);

        kk = 1:F_range;
        dlam_full = -2*sin(th_full*kk) * (kk(:) .* gco(2:F_range+1));

        eps_fold = interp1(th_e, eps_e, th_fold, 'spline');
        eps_full = sign(sin(th_full) + (sin(th_full)==0)) .* eps_fold;

        mu_full    = dlam_full .* eps_full;
        Fourier_mu = fourier_coefficients_spectral(mu_full, 2*F_range);

    % ==== Fix the sign with one test at moderate n ====

        %eig_f = widom_spectrum(a, n0, lambda_interval, Sidx, Cidx);
        eig_f = sort(eig(fourier_to_toeplitz(a, n0)));
        T_b   = fourier_to_toeplitz(FourierFP, n0);
        T_mu  = fourier_to_toeplitz(Fourier_mu, n0);

        d_plus  = sum(abs(eig_f(:) - sort(eig(T_b + T_mu/n0))));
        d_minus = sum(abs(eig_f(:) - sort(eig(T_b - T_mu/n0))));

        if d_minus < d_plus
            Fourier_mu = -Fourier_mu;
        end

        lam_e2 = lamS(th_e);
        c_inf  = abs(trapz(lam_e2, abs(eps_e)) / pi);

    % ==== Empirical cross-check of eps ====
        if opts.plot_phase_shift
            fs    = opts.fs;
            figure;
            hold on;
            plot(th_e, eps_e, 'b-', 'LineWidth', 2);
            grid on;
            box on;
            xlabel('$\theta$', 'Interpreter', 'latex', 'FontSize', fs);
            ylabel('$\varepsilon(\theta)$', 'Interpreter', 'latex', 'FontSize', fs);
            set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', fs);
            xlim([0 pi]);
            xticks(0:pi/4:pi);
            xticklabels({'$0$', '$\pi/4$', '$\pi/2$', '$3\pi/4$', '$\pi$'});
            set(gcf, 'Position', [100, 100, 500, 300]);
            hold off;
        end
end


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
% Determinants det(T_n(g) - lam) for all n in nvals, scaled by the first,
% from ONE banded no-pivot elimination (Sturm sequence): the k-th leading
% principal minor is the product of the first k pivots. Everything is
% tracked in log scale, so the output is O(1) as before.

    nmax = nvals(end);
    b    = min(numel(gco) - 1, nmax - 1);     % bandwidth

    c = zeros(nmax, 1);
    L = min(numel(gco), nmax);
    c(1:L) = gco(1:L);
    A = toeplitz(c);
    A(1:nmax+1:end) = A(1:nmax+1:end) - lam;  % T - lam*I

    % LAPACK-style clamp for (rare) near-singular leading minors
    pivmin = eps * (sum(abs(gco)) + abs(lam));

    ld = zeros(1, nmax);
    sg = zeros(1, nmax);
    lacc = 0;  sacc = 1;

    for k = 1:nmax
        d = A(k, k);
        if abs(d) < pivmin
            d = pivmin * (2*(d >= 0) - 1);
        end
        lacc  = lacc + log(abs(d));
        sacc  = sacc * (2*(d > 0) - 1);
        ld(k) = lacc;
        sg(k) = sacc;

        j2 = min(k + b, nmax);              % band-limited update, no fill-in
        if k < j2
            r = A(k, k+1:j2) / d;
            A(k+1:j2, k+1:j2) = A(k+1:j2, k+1:j2) - A(k+1:j2, k) * r;
        end
    end

    v = sg(nvals) .* exp(ld(nvals) - ld(nvals(1)));
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


function [Sidx, Cidx] = widom_subsets(m)
% All splittings of the 2m roots into m "inside" and m "outside".

    n2   = 2*m;
    Sidx = nchoosek(1:n2, m);
    nc   = size(Sidx, 1);

    % --- Complement of each row, vectorised ---
    mask = true(nc, n2);
    mask((Sidx - 1)*nc + (1:nc)') = false;     % clear the chosen entries
    Cidx = reshape(mod(find(mask.') - 1, n2) + 1, m, nc).';
end


function v = widom_det_seq(a, lambda, nvals, Sidx, Cidx)
% Scaled Widom determinants D_n / G^n for ALL n in nvals at once.
% Roots and combinatorial weights are computed once; n only enters
% through exp(n * logw), which becomes an outer product.

    m = (length(a)-1)/2;

    c = a;
    c(m+1) = c(m+1) - lambda;
    z = roots(c);
    [~, o] = sort(abs(z));
    z = z(o);

    lz   = log(z);
    lead = log((-1)^m * a(1));          % a(1) = a_{-m}

    zS  = z(Sidx);
    zC  = z(Cidx);
    slS = sum(lz(Sidx), 2);

    logw = lead + slS;
    D    = reshape(zS, [], m, 1) - reshape(zC, [], 1, m);
    lc   = m*slS - sum(sum(log(D), 3), 2);

    % --- n-dependent part: single outer product over all sizes ---
    w = logw - max(real(logw));         % column, length nchoosek(2m,m)
    v = real(sum(exp(lc + w * nvals(:).'), 1));   % 1 x numel(nvals)
end