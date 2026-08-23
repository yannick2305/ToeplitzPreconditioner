%{
    ----------------------------------------------------------------------
    Author(s):    [Yannick DE BRUIJN, Michael FLOATER, Erik HILTUNEN]
    Date:         [May 2026]
    Description:  [Numerical Preconditioner for Toeplitz matrices]
                  l1 spectral distance as a function of the matrix
                  dimension N, compared against the constant bound d_inf.
    ----------------------------------------------------------------------
%}
clc
clear;
close all;

% ==== Parameters ====
    m           = 5;            % Truncation size for a_k
    p           = 4.5;          % Decay rate upwards
    q           = 6.8;          % Decay rate downwards
    dim_range   = 5:60;         % Range of Toeplitz dimensions N to sweep
    num_lambda  = 30;           % Number of plotting points (50-100)
    fs          = 18;           % Fontsize for annotation
    normalize_l1 = false;       % true -> l1/N (mean spectral gap) instead of the raw sum

% ======================================================================
%  PART 1 : everything that is INDEPENDENT of the matrix dimension
%           (symbol, open limit, d_inf, GBZ, Fourier coefficients)
%           --> computed ONCE, before the dimension loop
% ======================================================================

% ==== Generate m-banded Dummy Toeplitz symbol ====
    col = zeros(m,1);
    row = zeros(1,m+1);

% --- Add noise to the coefficients ---
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

% ==== Approximate the open limit (dimension independent) ====
    [lambda_interval(1), lambda_interval(2)]         = open_limit(a);
    [lambda_interval_new(1), lambda_interval_new(2)] = open_limit_fast(a);

% ==== PRECOMPUTE d_inf ONCE ====
    d_inf = compute_cinf(a, lambda_interval(1), lambda_interval(2));

% --- compute_cinf returns a struct: extract the scalar bound c_inf ---
    if isstruct(d_inf)
        d_inf_val = d_inf.c_inf;
    else
        d_inf_val = d_inf;
    end
    d_inf_val = double(d_inf_val(1));   % force a plain scalar
    fprintf('Precomputed d_inf = %.8f\n', d_inf_val);

% ==== Approximate conjugate root set Lambda(f) ====
% --- Get sampling which is close to the DOS (Uniform spacing on Lambda(f)) ---
    t = linspace(0, 1, num_lambda);
    lambda_vals = lambda_interval(1) + (lambda_interval(2) - lambda_interval(1)) * (0.5 * (1 - cos(pi*t)));

% === Compute discrete points on GBZ ===
% --- Preallocate for all polynomial coefficients [P(z) = Q(z) - lambda*z^m] ---
    P_coeffs_matrix = repmat(a(:), 1, num_lambda);
% The coefficient of z^m is at position (m+1) in the array and subtract lambda
    P_coeffs_matrix(m + 1, :) = P_coeffs_matrix(m + 1, :) - lambda_vals;

% --- Compute all roots at once ---
    all_roots = zeros(2*m, num_lambda);
    for k = 1:num_lambda
        all_roots(:, k) = roots(P_coeffs_matrix(:, k));
    end

% --- Sort roots by magnitude for each lambda ---
    [~, sort_idx]    = sort(abs(all_roots), 1);
    linear_idx       = sort_idx + (0:num_lambda-1) * (2*m);
    all_roots_sorted = all_roots(linear_idx);

% --- Extract m-th and (m+1)-th roots ---
    candidate_roots = all_roots_sorted(m:m+1, :);

% --- Double check that the points lie on the GBZ ---
    tolerance = 1e-8;  % keep in mind rootsolver has its limits
    mod_n  = abs(candidate_roots(1, :));
    mod_n1 = abs(candidate_roots(2, :));
    similar_modulus = abs(mod_n - mod_n1) ./ max(mod_n, mod_n1) < tolerance;
    candidate_roots = candidate_roots(:, similar_modulus);

% --- Check if they have approximately the same modulus ---
    openLimit = [candidate_roots(1, :), candidate_roots(2, :)];
    phase     = angle(openLimit(:));
    [~, sortIdx] = sort(phase);
    openLimit_sorted = openLimit(sortIdx);

% --- Remove duplicates as they mess up the weights in the integral ---
    openLimit_sorted = merge_close_points(openLimit_sorted, 1e-1);

% ==================== For Michael ====================================
    nols  = length(openLimit_sorted);
    nolsh = nols/2;
    olsh  = openLimit_sorted(nolsh+2:nols);
    mnew  = num_lambda*5;
    mnew2 = mnew/2+1;
    interp_openLimit = resample2(olsh, mnew2);
% =====================================================================

% ==== Compute Fourier coefficients of f(p(z)) numerically ====
    k_values      = -m:m;
    phase_sorted  = angle(interp_openLimit(:));
    powers_matrix = interp_openLimit(:).^(-k_values);   % N x (2m+1) matrix

% --- Vectorized sum ---
    fp_values = powers_matrix * a(:);
    fp_values = real(fp_values);

% --- Compute the Fourier Transform of f(p(z)) ---
    F_range   = m + 15;   % truncate the Fourier coefficients to finite range
    FourierFP = fourier_coefficients_spectral(fp_values, F_range);

% ======================================================================
%  PART 2 : loop over the MATRIX DIMENSION and compute the l1 distance
% ======================================================================

    n_dims  = numel(dim_range);
    l1_vals = zeros(1, n_dims);

    for idx = 1:n_dims
        N = dim_range(idx);

        % --- Original Toeplitz matrix of dimension N ---
        T      = fourier_to_toeplitz(a, N);
        eigT   = sort(eig(T));

        % --- Toeplitz matrix for the deformed path, same dimension ---
        T_b    = fourier_to_toeplitz(FourierFP, N);
        eigT_b = sort(eig(T_b));

        % --- l1 (Wasserstein-1) distance between the two spectra ---
        l1_vals(idx) = sum(abs(eigT_b - eigT));

        if normalize_l1
            l1_vals(idx) = l1_vals(idx) / N;
        end

        fprintf('N = %3d   l1 = %.8e\n', N, l1_vals(idx));
    end

% ======================================================================
%  PART 3 : plot l1(N) together with the constant d_inf
% ======================================================================

    figure;
    hold on;
    plot(dim_range, l1_vals, 'b-o', 'LineWidth', 1.8, 'MarkerSize', 5, ...
         'MarkerFaceColor', 'w', 'DisplayName', '$\ell_1$ distance');
    plot(dim_range, d_inf_val*ones(size(dim_range)), 'r--', 'LineWidth', 2, ...
         'DisplayName', '$d_\infty$');
    grid on;
    box on;
    xlim([min(dim_range), max(dim_range)]);
    xlabel('Matrix dimension $n$', 'Interpreter', 'latex', 'FontSize', fs);
    %if normalize_l1
    %    ylabel('$\ell_1 / N$', 'Interpreter', 'latex', 'FontSize', fs);
    %else
    %    ylabel('$\ell_1(\sigma(\mathbf{T}_N), \sigma(\mathbf{T}^{b}_N))$', ...
    %          'Interpreter', 'latex', 'FontSize', fs);
   % end
    legend('Interpreter', 'latex', 'FontSize', fs, 'Location', 'best');
    set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', fs);
    set(gcf, 'Position', [100, 100, 500, 300]);
    hold off;


    %%

% --- Same plot on a logarithmic y-axis (often more informative) ---
    figure;
    semilogy(dim_range, l1_vals, 'b-o', 'LineWidth', 1.8, 'MarkerSize', 5, ...
             'MarkerFaceColor', 'w', 'DisplayName', '$\ell_1$ distance');
    hold on;
    semilogy(dim_range, d_inf_val*ones(size(dim_range)), 'r--', 'LineWidth', 2, ...
             'DisplayName', '$d_\infty$');
    grid on;
    box on;
    xlim([min(dim_range), max(dim_range)]);
    xlabel('Matrix dimension $N$', 'Interpreter', 'latex', 'FontSize', fs);
    ylabel('distance (log scale)', 'Interpreter', 'latex', 'FontSize', fs);
    legend('Interpreter', 'latex', 'FontSize', fs-4, 'Location', 'best');
    set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', fs);
    set(gcf, 'Position', [100, 100, 600, 380]);
    hold off;

% ======================================================================
%  PART 4 : convergence rate of l1(N) towards d_inf
%           expected model:  |l1(N) - d_inf| ~ C * N^(-1/3)
% ======================================================================

    Nv  = dim_range(:);
    err = abs(l1_vals(:) - d_inf_val);      % deviation from the limit
    sgn = sign(l1_vals(:) - d_inf_val);     % approach from above / below

    if mean(sgn) > 0
        fprintf('l1 approaches d_inf from ABOVE.\n');
    elseif mean(sgn) < 0
        fprintf('l1 approaches d_inf from BELOW.\n');
    else
        fprintf('l1 crosses d_inf inside the sampled range.\n');
    end

% --- Least-squares fit of the exponent on the asymptotic tail only ---
    fit_min_N = 20;                          % discard pre-asymptotic small N
    mask = (Nv >= fit_min_N) & (err > 0) & isfinite(err);

    if nnz(mask) < 3
        error('Not enough valid points for the rate fit - lower fit_min_N.');
    end

    Pfit  = polyfit(log(Nv(mask)), log(err(mask)), 1);
    alpha = Pfit(1);            % measured exponent
    Cfit  = exp(Pfit(2));       % measured prefactor

% --- Goodness of fit (R^2 in log-log coordinates) ---
    logerr_hat = polyval(Pfit, log(Nv(mask)));
    ss_res = sum((log(err(mask)) - logerr_hat).^2);
    ss_tot = sum((log(err(mask)) - mean(log(err(mask)))).^2);
    R2 = 1 - ss_res/ss_tot;

    fprintf('\n--- Convergence rate (fit over N >= %d) ---\n', fit_min_N);
    fprintf('  |l1 - d_inf| ~ %.6g * N^(%.4f)\n', Cfit, alpha);
    fprintf('  measured exponent : %.4f\n', alpha);
    fprintf('  expected exponent : %.4f   (-1/3)\n', -1/3);
    fprintf('  relative error    : %.2f %%\n', 100*abs(alpha-(-1/3))/(1/3));
    fprintf('  R^2 (log-log)     : %.4f\n', R2);

% --- Constant obtained by forcing the theoretical exponent -1/3 ---
    C_third = mean(err(mask) .* Nv(mask).^(1/3));
    fprintf('  C assuming N^(-1/3): %.6g\n\n', C_third);

% --- Local (pointwise) slopes: d log(err) / d log(N) ---
    valid = err > 0 & isfinite(err);
    Nvv   = Nv(valid);
    errv  = err(valid);
    local_slope = diff(log(errv)) ./ diff(log(Nvv));
    N_mid       = sqrt(Nvv(1:end-1) .* Nvv(2:end));   % geometric midpoints

% ---- Plot 1: log-log error with fitted and theoretical slopes ----
    figure;
    loglog(Nv, err, 'bo', 'MarkerSize', 6, 'LineWidth', 1.4, ...
           'MarkerFaceColor', 'w', 'DisplayName', '$|\ell_1(N)-d_\infty|$');
    hold on;
    Nplot = linspace(min(Nv), max(Nv), 200);
    loglog(Nplot, Cfit*Nplot.^alpha, 'b-', 'LineWidth', 1.8, ...
           'DisplayName', sprintf('fit: $N^{%.3f}$', alpha));
    loglog(Nplot, C_third*Nplot.^(-1/3), 'r--', 'LineWidth', 2, ...
           'DisplayName', '$N^{-1/3}$ reference');
    grid on; box on;
    xlabel('Matrix dimension $N$', 'Interpreter', 'latex', 'FontSize', fs);
    ylabel('$|\ell_1(N)-d_\infty|$', 'Interpreter', 'latex', 'FontSize', fs);
    legend('Interpreter', 'latex', 'FontSize', fs-4, 'Location', 'southwest');
    set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', fs);
    set(gcf, 'Position', [100, 100, 600, 380]);
    hold off;

% ---- Plot 2: local slope, should settle on -1/3 ----
    figure;
    semilogx(N_mid, local_slope, 'b-o', 'LineWidth', 1.6, 'MarkerSize', 5, ...
             'MarkerFaceColor', 'w', 'DisplayName', 'local slope');
    hold on;
    semilogx([min(N_mid) max(N_mid)], [-1/3 -1/3], 'r--', 'LineWidth', 2, ...
             'DisplayName', '$-1/3$');
    grid on; box on;
    xlabel('Matrix dimension $N$', 'Interpreter', 'latex', 'FontSize', fs);
    ylabel('$\mathrm{d}\log|\ell_1-d_\infty| / \mathrm{d}\log N$', ...
           'Interpreter', 'latex', 'FontSize', fs);
    ylim([-1.5, 0.5]);
    legend('Interpreter', 'latex', 'FontSize', fs-4, 'Location', 'best');
    set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', fs);
    set(gcf, 'Position', [100, 100, 600, 380]);
    hold off;

% ---- Plot 3: compensated error, should flatten to a constant ----
    figure;
    plot(Nv, err .* Nv.^(1/3), 'b-o', 'LineWidth', 1.6, 'MarkerSize', 5, ...
         'MarkerFaceColor', 'w', 'DisplayName', '$N^{1/3}|\ell_1-d_\infty|$');
    hold on;
    plot([min(Nv) max(Nv)], [C_third C_third], 'r--', 'LineWidth', 2, ...
         'DisplayName', sprintf('$C = %.4g$', C_third));
    grid on; box on;
    xlabel('Matrix dimension $N$', 'Interpreter', 'latex', 'FontSize', fs);
    ylabel('$N^{1/3}\,|\ell_1(N)-d_\infty|$', 'Interpreter', 'latex', 'FontSize', fs);
    legend('Interpreter', 'latex', 'FontSize', fs-4, 'Location', 'best');
    set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', fs);
    set(gcf, 'Position', [100, 100, 600, 380]);
    hold off;

%% --- Defining functions ---
function ck = fourier_coefficients_spectral(fp_values, K)
    M = length(fp_values/2);
% --- Use standard (uniform) FFT ---
    fft_result = fft(fp_values) / M;
% --- Extract coefficients ---
    ck = zeros(2*K+1, 1);
    ck(K+1)     = fft_result(1);
    ck(K+2:end) = fft_result(2:K+1);
    ck(1:K)     = fft_result(end-K+1:end);
end

function T = fourier_to_toeplitz(a, dimT)
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
    tol      = 1e-8;   % criterion: how flat is "zero"
    min_step = 1e-5;   % geometric: stop bisecting here (<< tol/|F'|)
    m = (length(a)-1)/2;
    coeff_Q = zeros(1, 2*m+1);
    for j = -m:m
        coeff_Q(m-j+1) = a(j+m+1);
    end
% Generate upper and lower bounds on the open limit
    f1 = sum(a);
    j = -m:m;
    f_minus1 = sum(a .* (-1).^j);
% --- Generate finite Toeplitz matrix ---
    T = fourier_to_toeplitz(a, 10);
    eigT = sort(eig(T));
    fmin = min(eigT);
    fmax = max(eigT);
% --- exact root computation ---
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
% inside:  confirmed inside the interval  (lambda0)
% outside: confirmed outside the interval (f1 or f_minus1)
        while abs(outside - inside) > min_step
            mid = (inside + outside) / 2;
            if in_interval(mid)
                inside = mid;
            else
                outside = mid;
            end
        end
        endpoint = inside;   % last confirmed inside point
    end
    lambdaR = find_endpoint(fmax, f1);
    lambdaL = find_endpoint(fmin, f_minus1);
end

function merged = merge_close_points(openLimit, tol)
    n = length(openLimit);
    if n == 0
        merged = [];
        return;
    end
    merged = openLimit;
% Check if first and last elements should merge
    if abs(openLimit(end) - openLimit(1)) <= tol
        merged(1) = mean([openLimit(1), openLimit(end)]);
        merged(end) = [];  % Remove last element
    end
% Check if two middle elements should merge
    if n >= 2
        mid1 = floor(n / 2);
        mid2 = mid1 + 1;
        if abs(merged(mid2) - merged(mid1)) <= tol
            merged(mid1) = mean([merged(mid1), merged(mid2)]);
            merged(mid2) = [];  % Remove second middle element
        end
    end
end

function [lambdaL, lambdaR] = open_limit_fast(a)
% OPEN_LIMIT_FAST  Endpoints of the (real) open/OBC limit set of a banded
% Toeplitz operator, computed directly instead of by bisection.
%
%   a = [a_{-m} ... a_0 ... a_m], symbol a(z) = sum_{j=-m}^{m} a_j z^j.
%
% Idea: lambda is an ENDPOINT of the limit set exactly when the polynomial
%   Q(z) = z^m (a(z) - lambda)
% has a double root among its middle pair, i.e. the m-th and (m+1)-th roots
% (sorted by modulus) coalesce. Writing out dQ/dz = 0 together with
% a(z0) = lambda gives simply a'(z0) = 0. So every endpoint is a critical
% value of the symbol: lambda = a(z0) with a'(z0) = 0. We compute all 2m
% critical points with ONE roots() call, keep the real critical values, and
% validate each with one roots() call (the collision must be the (m,m+1)
% pair by modulus, which is precisely the membership condition |r_m|=|r_{m+1}|).
%
% Cost: 1 + (#real critical values) calls to roots(), typically 3-6 total,
% versus ~2*log2(range/min_step) ~ 30+ calls for bisection - and the result
% is accurate to machine precision rather than min_step.
%
% Assumptions: a_m ~= 0 and a_{-m} ~= 0 (true bandwidth m); the limit set
% intersects the real axis and its real endpoints are what you want.
    a = a(:).';
    m = (numel(a) - 1) / 2;
    j = -m:m;
    P  = fliplr(a);        % descending coeffs of z^m * a(z)
    dP = fliplr(j .* a);   % descending coeffs of z^{m+1} * a'(z), degree 2m
    z0 = roots(dP);
    z0 = z0(abs(z0) > 1e-12 & isfinite(z0));   % discard spurious zero roots
    lambda = polyval(P, z0) ./ z0.^m;          % critical values a(z0)
% keep real candidates
    keep   = abs(imag(lambda)) < 1e-8 * max(1, abs(lambda));
    lambda = real(lambda(keep));
% validate: the double root must be the (m, m+1) modulus pair,
% i.e. lambda actually belongs to the limit set
    valid = false(size(lambda));
    for k = 1:numel(lambda)
        c = P;
        c(m+1) = c(m+1) - lambda(k);
        r = sort(abs(roots(c)));
        valid(k) = (r(m+1) - r(m)) <= 1e-6 * max(1, r(m+1));
    end
    lambda = lambda(valid);
    if isempty(lambda)
        error('open_limit_fast:noRealEndpoint', ...
              'No real endpoint found: the limit set may not meet the real axis.');
    end
    lambdaL = min(lambda);
    lambdaR = max(lambda);
end