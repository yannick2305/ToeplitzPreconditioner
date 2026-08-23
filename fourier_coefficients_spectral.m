
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