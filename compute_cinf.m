
function out = compute_cinf(cf)
% COMPUTE_CINF  Accurate computation of the asymptotic l1 spectral distance
%
%     c_inf = lim_{n->oo} d_sigma( T_n(f), T_n(f o p) )
%           = (1/pi) * int_0^pi |lambda'(a)| |DeltaPsi(a)| da
%           = (1/pi) * int_A^B  |DeltaPsi(lambda)| dlambda ,
%
% where DeltaPsi = Psi_g - Psi_f (mod pi) is the mismatch of the two
% quantisation phases,
%
%     Psi_h = (bandwidth)*alpha(lambda) - phi_h(lambda),
%     phi_h = (1/2) sum_{inner roots} arg( (zeta-z)/(conj(zeta)-z) )
%           + (1/2) sum_{outer roots} arg( (z-conj(zeta))/(z-zeta) ),
%
% built from the roots of z^bw*(h(z)-lambda), with zeta the upper root of
% the central conjugate pair.  Only small polynomial root-finding is
% needed: NO large-n eigenproblems, hence no spectral pollution and no
% variable precision -- double precision gives ~4 significant digits.
%
% Also returns the elementary lower-bound invariant
%     Delta2 = sum_k |k| ( g_k g_{-k} - a_k a_{-k} ),
% for which  c_inf >= |Delta2| / (2*Rhat).
%
% Reference values for m=7, kappa=3.5, rho=4.8 (350-bit certification):
%     c_inf  = 2.5993e-4       (measured d_sigma at n=200: 2.6007e-4)
%     Delta2 = +1.66545e-5

% ---------------- symbol (edit here) ----------------------------------
m    = floor(length(cf)/2);                        % bandwidth
%pdec = 4.5;                      % decay rate, a_k     = 1/(k+1)^pdec, k>0
%qdec = 4.8;                      % decay rate, a_{-k}  = 1/(k+1)^qdec
%kk   = 1:m;
%apos = 1./((kk+1).^pdec);
%aneg = 1./((kk+1).^qdec);
%a0   = 1;
% descending coefficients of z^m*(f(z)-lambda):  c(j) = a_{j-1-m}
%cf   = [aneg(end:-1:1), a0, apos];

% ---------------- numerical knobs -------------------------------------
NL   = 2000;                    % lambda samples for the band function
NF   = 2^13;                     % FFT length for the g_k
tg   = [28 34];                  % truncation bandwidths for g (stability pair)
NC   = 800;                      % lambda quadrature points for c_inf
edge = 1e-7;                     % relative margin excluded at band edges

% ---------------- 1. locate the band [A,B] ----------------------------
scan = linspace(0.3, 2.0, 4000);
inb  = false(size(scan));
for i = 1:numel(scan), inb(i) = in_band(scan(i), cf, m); end
ix = find(inb);
assert(~isempty(ix), 'no band found: widen the scan interval');
A = bisect(scan(ix(1)),  scan(max(ix(1)-1,1)),          cf, m);
B = bisect(scan(ix(end)),scan(min(ix(end)+1,numel(scan))), cf, m);
fprintf('band [A,B] = [%.8f, %.8f]\n', A, B);

% ---------------- 2. band function and Fourier coefficients g_k -------
lams = linspace(A+1e-12*(B-A), B-1e-12*(B-A), NL);
alph = zeros(1,NL);
for i = 1:NL
    z = central_pair(lams(i), cf, m);
    alph(i) = abs(angle(z));
end
[alph, ord] = sort(alph);  lam_s = lams(ord);
% append the exact edge values lambda(0), lambda(pi)
alph  = [0, alph, pi];
lam_s = [closest_edge(lam_s(1), A, B), lam_s, closest_edge(lam_s(end), A, B)];
[alph, iu] = unique(alph);  lam_s = lam_s(iu);
% clamped spline (lambda'(0)=lambda'(pi)=0, band function is even)
pp   = spline(alph, [0, lam_s, 0]);              % clamped ends
th   = 2*pi*(0:NF-1)/NF;
thf  = abs(mod(th+pi, 2*pi) - pi);               % fold to [0,pi]
gth  = ppval(pp, thf);
G    = fft(gth)/NF;
gk   = real(G(1:round(NF/2)));                   % g_k = g_{-k} real, gk(k+1)=g_k
fprintf('g_0 = %.10f  (a_0 = %g);   max imag/asym residue: %.1e\n', ...
        gk(1), a0, max(abs(imag(G(1:60)))));

% ---------------- 3. Delta_2 lower-bound invariant --------------------
K  = 120;
Delta2 = 2*sum((1:K).*gk(2:K+1).^2) - 2*sum(kk.*apos.*aneg);
Rhat   = max(sum(abs(cf)), gk(1)+2*sum(abs(gk(2:K+1))));
fprintf('Delta2 = %+.6e   =>   c_inf >= %.3e\n', Delta2, abs(Delta2)/(2*Rhat));

% ---------------- 4. phase mismatch and c_inf -------------------------
lamq = linspace(A+edge*(B-A), B-edge*(B-A), NC);
dPsi = zeros(numel(tg), NC);
for it = 1:numel(tg)
    t  = tg(it);
    cg = [gk(t+1:-1:2), gk(1), gk(2:t+1)];       % coeffs of z^t*(g_t - lam)
    for i = 1:NC
        Pf = psi_banded(lamq(i), cf, m);
        Pg = psi_banded(lamq(i), cg, t);
        dPsi(it,i) = mod(Pg - Pf + pi/2, pi) - pi/2;   % reduce mod pi
    end
end
stab = max(abs(dPsi(1,:) - dPsi(end,:)));
fprintf('t-stability of DeltaPsi (t=%d vs %d): %.2e\n', tg(1), tg(end), stab);
c_inf = trapz(lamq, abs(dPsi(end,:))) / pi;
fprintf('c_inf = %.6e\n', c_inf);

% ---------------- 5. optional validation at moderate n ----------------
n  = 40;                                          % double precision is safe here
Tf = toeplitz([a0, apos, zeros(1,n-m-1)], [a0, aneg, zeros(1,n-m-1)]);
kg = min(n-1, K);
Tg = toeplitz([gk(1), gk(2:kg+1), zeros(1,n-kg-1)]);
ds = sum(abs(sort(real(eig(Tf))) - sort(eig(Tg))));
fprintf('check: d_sigma at n=%d is %.4e  (should be near c_inf)\n', n, ds);

out = struct('c_inf',c_inf,'Delta2',Delta2,'A',A,'B',B,'gk',gk, ...
             'lam',lamq,'dPsi',dPsi(end,:),'stability',stab,'ds_check',ds);
end

% ======================================================================
function Psi = psi_banded(lam, c, bw)
% quantisation phase Psi = bw*alpha - phi for the symbol with descending
% coefficient vector c of z^bw*(h(z)-lambda-part); lam subtracted at centre.
c(bw+1) = c(bw+1) - lam;
r = roots(c);
[~, ix] = sort(abs(r));  r = r(ix);
z1 = r(bw);  z2 = r(bw+1);
if imag(z1) > 0, zeta = z1; else, zeta = z2; end
inner = r(1:bw-1);  outer = r(bw+2:end);
alpha = angle(zeta);
phi   = 0.5*( sum(angle((zeta - inner)./(conj(zeta) - inner))) ...
            + sum(angle((outer - conj(zeta))./(outer - zeta))) );
Psi = bw*alpha - phi;
end

function z = central_pair(lam, cf, m)
c = cf;  c(m+1) = c(m+1) - lam;
r = roots(c);
[~, ix] = sort(abs(r));  r = r(ix);
if imag(r(m)) > 0, z = r(m); else, z = r(m+1); end
end

function tf = in_band(lam, cf, m)
c = cf;  c(m+1) = c(m+1) - lam;
r = roots(c);
[~, ix] = sort(abs(r));  r = r(ix);
z1 = r(m); z2 = r(m+1);
tf = abs(imag(z1)) > 1e-12*abs(z1) && ...
     abs(abs(z1) - abs(z2)) < 1e-9*abs(z1);      % conjugate pair
end

function x = bisect(xin, xout, cf, m)
for it = 1:80
    xm = (xin + xout)/2;
    if in_band(xm, cf, m), xin = xm; else, xout = xm; end
end
x = xin;
end

function y = closest_edge(y0, A, B)
if abs(y0 - A) < abs(y0 - B), y = A; else, y = B; end
end
