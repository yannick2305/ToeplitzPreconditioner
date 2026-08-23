%{
    ----------------------------------------------------------------------
    Author(s):    [Yannick DE BRUIJN, Michael FLOATER, Erik HILTUNEN]
    Date:         [August 2026]
    Description:  [Run File]
    ----------------------------------------------------------------------
%}

clear;
close all;
clc

% ==== Parameters ====
    m = 3;      % Truncation size for a_k
    p = 3;    % Decay rate upwards
    q = 4;    % Decay rate downwards
    DimT = 80;  % matrix dimension

    Show_the_GBZ     = true;
    Show_phase_shift = true;

    a = generate_symbol(m, p, q);
    Ta = fourier_to_toeplitz(a, DimT);
    eigTa = sort(eig(Ta));

% === First order expansion and GBZ ===
    F1 = first_order_expansion(a, Show_the_GBZ);
    T1 = fourier_to_toeplitz(F1, DimT);
    eigT1 = sort(eig(T1));

% === Second order expansion ===
    [F_1, F_2, c_inf] = first_and_second_order_expansion(a, plot_phase_shift=true, n_eps=200);
    T_1 = fourier_to_toeplitz(F_1, DimT);
    T_2 = fourier_to_toeplitz(F_2, DimT);
    T_Sec = T_1 + T_2/DimT;
    eigT_Sec = sort(eig(T_Sec));

    % --- Compare the spectra ---
    l1 = sum(abs(eigTa - eigT1));
    fprintf('Waaserstein First order : %f\n', l1);

    l2 = sum(abs(eigTa - eigT_Sec));
    fprintf('Waaserstein Second order: %f\n', l2);


% === Convergence of the Wasserstein Distance ===

    % --- Sweep over matrix sizes ---
    DimList = unique(round(logspace(log10(10), log10(100), 12)));
    d_unc   = zeros(size(DimList));
    d_cor   = zeros(size(DimList));

    for k = 1:numel(DimList)
        n = DimList(k);

        eig_f = sort(eig(fourier_to_toeplitz(a, n)));   

        % first order only
        eigT1    = sort(eig(fourier_to_toeplitz(F1, n)));
        d_unc(k) = sum(abs(eig_f - eigT1));

        % first + second order
        T_1      = fourier_to_toeplitz(F_1, n);
        T_2      = fourier_to_toeplitz(F_2, n);
        eigT_Sec = sort(eig(T_1 + T_2/n));
        d_cor(k) = sum(abs(eig_f - eigT_Sec));
    end

    % --- Plot the Convergence ---
    fs = 16;
    figure;
    loglog(DimList, d_unc, 'k-o', 'LineWidth', 1.5);
    hold on;
    loglog(DimList, d_cor,                         'b-s', 'LineWidth', 1.5);
    loglog(DimList, abs(c_inf)*ones(size(DimList)), 'g:', 'LineWidth', 1.5);
    loglog(DimList, d_cor(1)*DimList(1)./DimList,  'r--', 'LineWidth', 1.5);
    grid on;
    box on;
    xlabel('Matrix size $n$', 'Interpreter', 'latex', 'FontSize', fs+2);
    ylabel('$d_\sigma$', 'Interpreter', 'latex', 'FontSize', fs+2);
    legend({'$d_\sigma(\mathbf{T}_n(f), \mathbf{T}_n(f\circ p))$', ...
            '$d_\sigma(\mathbf{T}_n(f), \mathbf{T}_n(f\circ p) + \frac{1}{n}\mathbf{T}_n(\lambda''\varepsilon))$', ...
            '$c_\infty$ prediction', ...
            '$\mathcal{O}(n^{-1})$'}, ...
            'Interpreter', 'latex', 'Location', 'southwest');
    set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', fs);
    ylim([d_cor(end)*0.4, c_inf*1.3]);
    set(gcf, 'Position', [100, 100, 500, 300]);
    hold off;


  



