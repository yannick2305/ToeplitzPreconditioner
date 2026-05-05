%{
    ----------------------------------------------------------------------
    Author(s):    [Erik Orvehed HILTUNEN , Yannick DE BRUIJN]
    Date:         [December 2025]
    Description:  [Asymptotic Similarity Transform for Toeplitz matrices]
    ----------------------------------------------------------------------
%}

clc
clear;
close all;


% ==== Parameters ====
    m = 9;              % Truncation size for a_k
    p = 3.5;            % Decay rate upwards
    q = 4.8;            % Decay rate downwards
    DimT = 600;         % Dimension of finite Toeplitz matrix to simulate open limit
    num_lambda = 40;    % Number of plotting points (50-300)
    fs = 18;            % Fontsize for annotation
    

% ==== Generate m-bsned Dummy Toeplitz Matrix ====
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
       
    % --- Generate finite Toeplitz matrix ---
    T = fourier_to_toeplitz(a, DimT);
    eigT = sort(eig(T));


% ==== Approximate the open limit ====
    % --- Get initial guess inside the open limit ---
    lambda_start = ( min(real(eigT)) + max(real(eigT)) ) / 2;

    % --- Adaptive Computation of the Open Limit ---
     [lambda_interval(1), lambda_interval(2)] = open_limit(a);
     %[lambda_interval_1(1), lambda_interval_1(2)] = open_limit_warm(a);


% ==== Approximate conjugate root set Λ(f) ====
    % --- Get sampling which is close to the DOS (Uniform spacing on Λ(f) ---
    t = linspace(0, 1, num_lambda);
    lambda_vals = lambda_interval(1) + (lambda_interval(2) - lambda_interval(1)) * (0.5 * (1 - cos(pi*t)));

    % --- Plot the Density of states of the sampling ---
    %{
    histogram(lambda_vals, floor(num_lambda/40), 'Normalization', 'pdf');
    xlabel('Eigenvalue');
    ylabel('Density of States');
    %}


% === Compute discrete points on GBZ ===
    % --- Preallocate for all polynomial coefficients [P(z) = Q(z) - lambda*z^m] ---
    P_coeffs_matrix = repmat(a(:), 1, num_lambda);
    
    % The coefficient of z^m is at position (m+1) in the array and substract lambda
    P_coeffs_matrix(m + 1, :) = P_coeffs_matrix(m + 1, :) - lambda_vals;
    
    % --- Compute all roots at once ---
    all_roots = zeros(2*m, num_lambda);
    for k = 1:num_lambda
        all_roots(:, k) = roots(P_coeffs_matrix(:, k));
    end
    
    % --- Sort roots by magnitude for each lambda ---
    [~, sort_idx] = sort(abs(all_roots), 1);
    linear_idx = sort_idx + (0:num_lambda-1) * (2*m);
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

    phase = angle(openLimit(:));
    [~, sortIdx] = sort(phase);
    openLimit_sorted = openLimit(sortIdx);

    % --- Remove duplicated as they mess up the weights in integral ---
    openLimit_sorted = merge_close_points(openLimit_sorted, 1e-4);

    
% === Spline to get continuous GBZ ===
    chord = resampleCurveByAngle(openLimit_sorted, 1000);

    % --- Plot the set Λ(f) ---
    %
        wraparound_OpenLimit = [openLimit_sorted, openLimit_sorted(1)];

        figure;
        % --- discrete GBZ ---
        plot(real(wraparound_OpenLimit), imag(wraparound_OpenLimit), 'x-', 'LineWidth', 2.5)
        hold on;
        % --- Spline Interpolated continous GBZ ---
        plot(real(chord), imag(chord), 'm-', 'LineWidth', 2.5)
        plot(real(chord), imag(chord), 'mx', 'LineWidth', 2.5)
        
        % ---- unit circle ----
        theta = linspace(0, 2*pi, 300);
        plot(cos(theta), sin(theta), 'r-', 'LineWidth', 2);
        %legend({'$\Lambda(f)$', 'Unit torus', ''}, 'Interpreter', 'latex', 'FontSize', 18);
        set(gcf, 'Position', [100, 100, 300, 300]);
        set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 18);
        xlabel('$\mathrm{Re}$', 'Interpreter', 'latex', 'FontSize', 18);
        ylabel('$\mathrm{Im}$', 'Interpreter', 'latex', 'FontSize', 18); 
        grid on;
        axis equal;
        box on;
        hold off;
    %}


% ==== Compute Fourier coefficients of f(p(z)) numerically  ====
    % --- Evaluate f(p(z)) on the torus ---
    k_values = -m:m;
    openLimit_sorted = chord;
    phase_sorted = angle(openLimit_sorted(:));
    powers_matrix = openLimit_sorted(:).^(-k_values);  % N x (2n+1) matrix
    
    % --- Vectorized sum ---
    fp_values = powers_matrix * a(:);

    % --- Clean up data ---
    fp_values = real(fp_values);

    % --- Wrap around ---
    phase_sorted = [-pi; phase_sorted; pi];
    fp_values    = [ fp_values(end); fp_values; fp_values(end)];

    % --- Plot the function f(p(torus))
    %{
    figure;
    plot(phase_sorted, fp_values);
    hold off;
    %}

    % --- Compute the Fourier Transform of f(p(z)) ---
    F_range = m+18;
    FourierFP = fourier_coefficients_spectral(fp_values, F_range);

    % --- Plot the decay in the Fourier Coefficients ---
    %{
    figure;
    plot(1:length(FourierFP), log(abs(FourierFP)));
    %}

% ==== Quasi Similarity Transformed Toeplitz matrix ====
    % --- Toeplitz matrix for deformed path ---
    T_b = fourier_to_toeplitz(FourierFP, DimT);
    eigT_b = sort(eig(T_b));

    % --- Plot the eigenvalues before and after asymptotic Similarity transform ---
    %
    figure;
    plot(real(eigT), imag(eigT), 'bx', 'MarkerSize', 8, 'LineWidth', 1.5);
    hold on;
    plot(real(eigT_b), imag(eigT_b), 'ro', 'MarkerSize', 8, 'LineWidth', 1.5);
    grid on;
    box on;
    xlim([0.96*min(real(eigT)), 1.03*max(real(eigT))])
    ylim([-0.005, 0.005])
    xlabel('$\mathrm{Re}(\sigma(\mathbf{T}_N))$', 'Interpreter', 'latex', 'FontSize', 14);
    ylabel('$\mathrm{Im}(\sigma(\mathbf{T}_N))$', 'Interpreter', 'latex', 'FontSize', 14); 
    set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 18);
    set(gcf, 'Position', [100, 100, 500, 300]); 
    axis equal;
    hold off;
    %}

% ==== Hirschman density of states ====
    %{
    lambda_interval = linspace(lambda_interval(1)-0.1, lambda_interval(2)+0.1, 100);
    deriv_sums = zeros(size(lambda_interval));
    
    for i = 1:length(lambda_interval)
        deriv_sums(i) = 1/(2*pi) * 1/g(lambda_interval(i), a) * compute_g_derivative_sum(lambda_interval(i), a);
        deriv_sums(i) = min(20, deriv_sums(i));
    end
    %}

    % --- Plot the DOS against the empirical measure ---
    %{
    figure;
    histogram(eigT_b, floor(DimT/3), 'Normalization', 'pdf');
    hold on;
    plot(lambda_interval, deriv_sums, 'LineWidth', 2);
    set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 18);
    set(gcf, 'Position', [100, 100, 600, 250]); 
    xlabel('Eigenvalue', 'Interpreter', 'latex', 'FontSize', 14);
    ylabel('Density of States', 'Interpreter', 'latex', 'FontSize', 14); 
    box on;
    hold off;
    %}


%% --- Defining functions ---


function ck = fourier_coefficients_spectral(fp_values, K)

    M = length(fp_values/2);


    fp_uniform = fp_values;
    
    % --- Use standard (uniform) FFT ---
    fft_result = fft(fp_uniform) / M;
    
    % --- Extract coefficients ---
    ck = zeros(2*K+1, 1);
    ck(K+1) = fft_result(1);
    ck(K+2:end) = fft_result(2:K+1);
    ck(1:K) = fft_result(end-K+1:end);
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
    min_step = 1e-5;   % geometric: stop bisecting here (≪ tol/|F'|)

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

    function F = diff_mod(lambda)
        c      = coeff_Q;
        c(m+1) = c(m+1) - lambda;
    
        A = compan(c);               % build companion matrix directly
    
        opts.tol = 1e-6;             % loose tolerance
        opts.p   = min(3*(m+1), 2*m); % Krylov subspace size
        d = eigs(A, m+1, 'sm', opts); % only the m+1 smallest-magnitude eigenvalues
    
        mags = sort(abs(d));
        F    = mags(m+1) - mags(m);
    end

    %{
    % --- exact root computation ---
    function F = diff_mod(lambda)
        c = coeff_Q;
        c(m+1) = c(m+1) - lambda;
        r = roots(c);
        r = sort(abs(r));
        F = abs(r(m+1)) - abs(r(m)) - abs(imag(r(m+1)));
    end
    %}

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


%{
function [lambdaL, lambdaR] = open_limit_warm(a)

    tol      = 1e-8;
    min_step = 1e-10;
    m        = (length(a)-1)/2;

    coeff_Q = zeros(1, 2*m+1);
    for j = -m:m
        coeff_Q(m-j+1) = a(j+m+1);
    end

    f1       = sum(a);
    j        = -m:m;
    f_minus1 = sum(a .* (-1).^j);

    % --- Generate finite Toeplitz matrix ---
    T    = fourier_to_toeplitz(a, 10);
    eigT = sort(eig(T));
    fmin = min(eigT);
    fmax = max(eigT);

    % eigs options — built once, reused every call
    opts.tol = 1e-6;
    opts.p   = min(3*(m+1), 2*m);

    v0_cache = [];

    function F = diff_mod(lambda)
        c      = coeff_Q;
        c(m+1) = c(m+1) - lambda;
        A      = compan(c);

        if ~isempty(v0_cache)
            opts.v0 = v0_cache;
        end

        [V, D]   = eigs(A, m+1, 'sm', opts);
        v0_cache = V(:, 1);            % best Ritz vector for next call

        mags = sort(abs(diag(D)));
        F    = mags(m+1) - mags(m);
    end

    function inside = in_interval(lambda)
        inside = abs(diff_mod(lambda)) <= tol;
    end

    function endpoint = find_endpoint(inside_pt, outside_pt)
        while abs(outside_pt - inside_pt) > min_step
            mid = (inside_pt + outside_pt) / 2;
            if in_interval(mid)
                inside_pt = mid;
            else
                outside_pt = mid;
            end
        end
        endpoint = inside_pt;
    end

    % Prime cache at fmax, search right
    v0_cache = [];
    diff_mod(fmax);
    lambdaR = find_endpoint(fmax, f1);

    % Reset cache, prime at fmin, search left
    v0_cache = [];
    diff_mod(fmin);
    lambdaL = find_endpoint(fmin, f_minus1);

end
%}



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


function g_val = g(lambda, a, coeff_base)
% USAGE:
%   g_val = g(lambda, a)              % standalone
%   g_val = g(lambda, a, coeff_base)  % fast path: pass precomputed base

    m = (length(a) - 1) / 2;

    % Precompute only if not supplied (avoids recomputation in hot loops)
    if nargin < 3
        coeff_base = fliplr(a(:).');   % replaces the for loop
    end

    % Shift constant term for this lambda
    c      = coeff_base;
    c(m+1) = c(m+1) - lambda;

    % Partial sort: only need the m largest magnitudes
    r_abs        = abs(roots(c));
    larger_m_abs = maxk(r_abs, m);

    g_val = abs(a(end)) * prod(larger_m_abs);
end


function deriv_sum = compute_g_derivative_sum(lambda_real, a, h)

    if nargin < 3
        h = 1e-5;
    end

    lambda = complex(lambda_real, 0);

    % Evaluate once per shift
    g_plus   = g(lambda + 1i*h, a);
    g_center = g(lambda, a);
    g_minus  = g(lambda - 1i*h, a);

    deriv_sum = (g_plus - 2*g_center + g_minus) / h;
end



function zUniform = resampleCurveByAngle(z, N)
    % Interpolates using chord-length parameterization for smoothness,
    % but resamples so that angle(zUniform) is uniformly spaced.
    % Assumes the curve is star-shaped w.r.t. the origin.
    
    x = real(z(:));
    y = imag(z(:));
    n = numel(x);
    
    % 1. Chord-length parameterization
    dx = diff([x; x(1)]);
    dy = diff([y; y(1)]);
    ds = sqrt(dx.^2 + dy.^2);
    t  = [0; cumsum(ds)];
    t  = t(1:end-1) / t(end);   % n values in [0,1)
    
    % 2. Periodic extension for smooth spline
    nExt = 3;
    tExt = [ t(end-nExt+1:end)-1 ; t ; t(1:nExt)+1 ];
    xExt = [ x(end-nExt+1:end)   ; x ; x(1:nExt)   ];
    yExt = [ y(end-nExt+1:end)   ; y ; y(1:nExt)   ];
    
    % 3. Evaluate on a dense grid to get a finely sampled curve
    M      = max(10000, 100*n);
    tDense = linspace(0, 1, M+1)';
    xD     = interp1(tExt, xExt, tDense, 'spline');
    yD     = interp1(tExt, yExt, tDense, 'spline');
    zDense = xD + 1i*yD;
    
    % 4. Compute the unwrapped angle along the dense curve
    thetaDense = unwrap(angle(zDense));
    
    % Ensure monotonicity (should be guaranteed for a star-shaped curve,
    % but fix any tiny numerical non-monotone steps just in case)
    thetaDense = cummax_monotone(thetaDense);
    
    % 5. Define uniform target angles over exactly one full revolution
    theta0    = thetaDense(1);
    thetaUnif = linspace(theta0, theta0 + 2*pi, N+1)';
    thetaUnif = thetaUnif(1:end-1);   % drop duplicate endpoint
    
    % 6. Invert: find the chord-length parameter t where angle == thetaUnif
    tUnif = interp1(thetaDense, tDense, thetaUnif, 'pchip');
    
    % 7. Evaluate the spline at those t values
    xU       = interp1(tExt, xExt, tUnif, 'spline');
    yU       = interp1(tExt, yExt, tUnif, 'spline');
    zUniform = xU + 1i*yU;
end

% Helper: force strict monotone increase (fixes tiny numerical glitches)
function y = cummax_monotone(x)
    y = x;
    for k = 2:numel(x)
        if y(k) <= y(k-1)
            y(k) = y(k-1) + 1e-14;
        end
    end
end
