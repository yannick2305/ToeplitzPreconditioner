%{
    ----------------------------------------------------------------------
    Author(s):    [Yannick DE BRUIJN, Michael FLOATER, Erik HILTUNEN]
    Date:         [May 2026]
    Description:  [Numerical Preconditioner for Toeplitz matrices]
    ----------------------------------------------------------------------
%}

clc
clear;
close all;


% ==== Parameters ====
    m = 7;              % Truncation size for a_k
    p = 3.5;            % Decay rate upwards
    q = 4.8;            % Decay rate downwards
    DimT = 40;          % Dimension of finite Toeplitz matrix to simulate open limit
    num_lambda = 30;    % Number of plotting points (50-100)
    fs = 18;            % Fontsize for annotation
    

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
       
    % --- Generate finite Toeplitz matrix ---
    T = fourier_to_toeplitz(a, DimT);
    eigT = sort(eig(T));


% ==== Approximate the open limit ====
     [lambda_interval(1), lambda_interval(2)] = open_limit(a);


% ==== Approximate conjugate root set Λ(f) ====
    % --- Get sampling which is close to the DOS (Uniform spacing on Λ(f)) ---
    t = linspace(0, 1, num_lambda);
    lambda_vals = lambda_interval(1) + (lambda_interval(2) - lambda_interval(1)) * (0.5 * (1 - cos(pi*t)));


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
    openLimit_sorted = merge_close_points(openLimit_sorted, 1e-1);


    % ==================== For Michael ====================================

    % openLimit_sorted are the complex values, sorted in their non-uniform  
    % phase. These points are to be interpolated and resampled uniformly to
    % pass them into an FFT.

    nols  = length(openLimit_sorted);
    nolsh = nols/2;
    olsh  = openLimit_sorted(nolsh+2:nols);
    mnew  = num_lambda*5;
    mnew2 = mnew/2+1;

 
    interp_openLimit = resample2(olsh, mnew2);
    length(interp_openLimit)

    % =====================================================================

  
    % --- Plot the set Λ(f) ---
    %
        wraparound_Int_GBZ = [interp_openLimit, interp_openLimit(1)];

        figure;
        % --- discrete GBZ ---
        plot(real(openLimit_sorted), imag(openLimit_sorted), 'bx', 'LineWidth', 2.5)
        hold on;
        % --- Interpolated continous GBZ ---
        plot(real(wraparound_Int_GBZ), imag(wraparound_Int_GBZ), 'm-', 'LineWidth', 2.5)

        % ---- unit circle ----
        theta = linspace(0, 2*pi, 300);
        plot(cos(theta), sin(theta), 'r-', 'LineWidth', 2);

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

    % --- Evaluate f(p(z)) on the torus, where p is interpolated ---
    k_values = -m:m;
    phase_sorted = angle(interp_openLimit(:));
    powers_matrix = interp_openLimit(:).^(-k_values);  % N x (2n+1) matrix
    
    % --- Vectorized sum ---
    fp_values = powers_matrix * a(:);

    % --- Clean up data ---
    fp_values = real(fp_values);

    % --- Compute the Fourier Transform of f(p(z)) ---
    F_range = m + 15; % truncate the Fourier coefficients to finite range
    FourierFP = fourier_coefficients_spectral(fp_values, F_range);

    % --- Plot the decay in the Fourier Coefficients ---
    %{
    figure;
    plot(1:length(FourierFP), log(abs(FourierFP)));
    %}    

    % --- Toeplitz matrix for deformed path ---
    T_b = fourier_to_toeplitz(FourierFP, DimT);
    eigT_b = sort(eig(T_b));

    
    l1 = sum(abs(eigT_b - eigT));
    fprintf('Waaserstein distance between spectra: %f\n', l1);

    
    % --- Plot the eigenvalues before and after asymptotic Similarity transform ---
    %{
    figure;
    hold on;
    plot(real(eigT),   imag(eigT),   'kx', 'MarkerSize', 8, 'LineWidth', 1.5);
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


