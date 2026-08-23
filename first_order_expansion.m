function FourierFP = first_order_expansion(a, showGBZ)


    %{
        ----------------------------------------------------------------------
        Author(s):    [Yannick DE BRUIJN, Michael FLOATER, Erik HILTUNEN]
        Date:         [August 2026]
        Description:  [First order expansion]
        ----------------------------------------------------------------------
    %}

    arguments
        a
        showGBZ logical = false   % optional, defaults to false
    end

    num_lambda = 40;    % Number of plotting points (50-100)
    m = floor(length(a)/2);

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

    % ==== Compute Fourier coefficients of f(p(z)) numerically  ====

    % --- Evaluate f(p(z)) on the torus, where p is interpolated ---
    k_values = -m:m;
   
    powers_matrix = interp_openLimit(:).^(-k_values);  % N x (2n+1) matrix
    
    % --- Vectorized sum ---
    fp_values = powers_matrix * a(:);

    % --- Clean up data ---
    fp_values = real(fp_values);

    % --- Compute the Fourier Transform of f(p(z)) ---
    F_range = m + 15; % truncate the Fourier coefficients to finite range
    FourierFP = fourier_coefficients_spectral(fp_values, F_range);


    % --- Plot the generalised Brillouin zone ---
    if showGBZ
        
        wraparound_Int_GBZ = [interp_openLimit, interp_openLimit(1)];

        figure;
        % --- discrete GBZ ---
        %plot(real(openLimit_sorted), imag(openLimit_sorted), 'bx', 'LineWidth', 2.5)
        hold on;
        % --- Interpolated continous GBZ ---
        hGBZ = plot(real(wraparound_Int_GBZ), imag(wraparound_Int_GBZ), 'r-', 'LineWidth', 2.5);
        % ---- unit circle ----
        theta = linspace(0, 2*pi, 300);
        plot(cos(theta), sin(theta), 'k--', 'LineWidth', 2);
        legend(hGBZ, {'GBZ'}, ...
            'Interpreter', 'latex', 'FontSize', 14, 'Location', 'northeast');
        set(gcf, 'Position', [100, 100, 300, 300]);
        set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 18);
        xlabel('$\mathrm{Re}$', 'Interpreter', 'latex', 'FontSize', 18);
        ylabel('$\mathrm{Im}$', 'Interpreter', 'latex', 'FontSize', 18);
        grid on;
        axis equal;
        box on;
        ylim([min(imag(wraparound_Int_GBZ))*1.3, max(imag(wraparound_Int_GBZ))*1.3]);
        hold off;

    end

end