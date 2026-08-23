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