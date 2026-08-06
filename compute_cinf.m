function out = compute_cinf(cf, A, B, opts)
% COMPUTE_CINF  Asymptotic l1 spectral distance c_inf = lim_n d_sigma.
%
%   out = compute_cinf(cf)            % locates the band itself
%   out = compute_cinf(cf, A, B)      % band [A,B] supplied (faster)
%   out = compute_cinf(cf, A, B, opts)
%
%   cf   descending coefficients (a_{-m},...,a_0,...,a_m) of z^m*f_m(z)
%   A,B  endpoints of the open limit (e.g. from the bisection algorithm)
%   opts struct, optional fields:
%          .t      truncation bandwidth of g = f o p        (default 30)
%          .NL     band-function samples                    (default 100)
%          .NF     FFT length for the g_k                   (default 2^12)
%          .ng     Gauss nodes per smooth panel             (default 12)
%          .npan   panels per smooth piece                  (default 2)
%          .check  run the n=40 eigenvalue cross-check      (default true)
%          .verbose                                          (default true)
%
%   Returns a struct with c_inf, Delta2, the bracket, the phase profile,
%   and an internal error estimate.
%
% METHOD.  c_inf = (1/pi) int_A^B |DeltaPsi(lambda)| dlambda, where
% DeltaPsi = Psi_g - Psi_f is the mismatch of the two quantisation phases.
% Three devices make this cheap and accurate:
%
%  (1) CANCELLATION-FREE PHASES.  Writing zeta for the upper central root
%      and zb = conj(zeta), each phase is reorganised as
%
%        Psi_h = alpha + (1/2) sum_inner arg( zeta*(zb-z) / (zb*(zeta-z)) )
%                      - (1/2) sum_outer arg( (w-zb)/(w-zeta) ),
%
%      in which every summand tends to 0 as an inner root tends to 0 or an
%      outer root to infinity.  The leading alpha cancels in DeltaPsi, so
%      no large quantities are subtracted, the result is independent of the
%      bandwidth by construction, and NO mod-pi reduction is needed.
%
%  (2) PALINDROMIC REDUCTION (the main speed-up).  Since g_k = g_{-k}, the
%      polynomial z^t*(g_t(z)-lambda) is palindromic: its roots come in
%      pairs (z,1/z), and with u = z + 1/z it collapses to a degree-t
%      polynomial whose CHEBYSHEV coefficients are exactly (g_0-lambda,
%      2*g_1, ..., 2*g_t).  Its roots are the eigenvalues of a t-by-t
%      colleague matrix -- an 8x saving over a degree-2t root solve, and
%      better conditioned (no basis conversion).
%
%  (3) QUADRATURE.  lambda = (A+B)/2 - (B-A)/2*cos(theta) absorbs the
%      square-root behaviour of DeltaPsi at the band edges exactly; the
%      finitely many sign changes of DeltaPsi (the kinks of |DeltaPsi|) are
%      located by bisection and used as panel boundaries.  Composite
%      Gauss-Legendre on the smooth panels then converges spectrally:
%      ~10^-12 with about 120 integrand evaluations.

if nargin < 4, opts = struct(); end
if nargin < 3, A = []; B = []; end
def = @(f,v) (isfield(opts,f) && ~isempty(opts.(f)))*1;
t       = 30;    if def('t',1),       t       = opts.t;       end
NL      = 100;   if def('NL',1),      NL      = opts.NL;      end
NF      = 2^12;  if def('NF',1),      NF      = opts.NF;      end
ng      = 12;    if def('ng',1),      ng      = opts.ng;      end
npan    = 2;     if def('npan',1),    npan    = opts.npan;    end
docheck = true;  if def('check',1),   docheck = opts.check;   end
verbose = true;  if def('verbose',1), verbose = opts.verbose; end

cf = cf(:).';                      % row, descending
m  = floor(numel(cf)/2);
a0 = cf(m+1);  apos = cf(m+2:2*m+1);  aneg = cf(m:-1:1);

% ---------- 0. band -----------------------------------------------------
if isempty(A) || isempty(B)
    [A, B] = locate_band(cf, m);
end
if verbose, fprintf('band [A,B] = [%.10f, %.10f]\n', A, B); end
mid = (A+B)/2;  hw = (B-A)/2;
lam_of_th = @(th) mid - hw*cos(th);
th_of_lam = @(l)  acos(max(-1, min(1, (mid-l)/hw)));

% ---------- 1. band function and Fourier coefficients g_k ---------------
% Chebyshev-spaced lambda nodes => nearly uniform alpha nodes, because
% lambda - edge ~ alpha^2 at both ends.  This is what lets NL stay small.
thL  = linspace(0, pi, NL);
lams = lam_of_th(thL);
lams([1 end]) = [A + 1e-13*(B-A), B - 1e-13*(B-A)];
alph = zeros(1, NL);
for i = 1:NL
    alph(i) = abs(angle(central_root(lams(i), cf, m)));
end
[alph, ord] = sort(alph);  lam_s = lams(ord);
alph  = [0, alph, pi];
lam_s = [edge_of(lam_s(1),A,B), lam_s, edge_of(lam_s(end),A,B)];
[alph, iu] = unique(alph);  lam_s = lam_s(iu);
pp  = spline(alph, [0, lam_s, 0]);                  % clamped: lambda'=0 at ends
th  = 2*pi*(0:NF-1)/NF;
gth = ppval(pp, abs(mod(th+pi, 2*pi) - pi));        % even extension
G   = fft(gth)/NF;
gk  = real(G(1:round(NF/2)));                       % gk(k+1) = g_k = g_{-k}

% ---------- 2. elementary invariants ------------------------------------
K      = min(120, numel(gk)-1);
Delta2 = 2*sum((1:K).*gk(2:K+1).^2) - 2*sum((1:m).*apos.*aneg);
Rhat   = max(sum(abs(cf)), gk(1) + 2*sum(abs(gk(2:K+1))));
floorb = abs(Delta2)/(2*Rhat);
if verbose
    fprintf('g_0 = %.10f (a_0 = %g),  |g_%d| = %.1e\n', gk(1), a0, t, abs(gk(t+1)));
    fprintf('Delta2 = %+.6e  =>  c_inf >= %.3e\n', Delta2, floorb);
end

% ---------- 3. locate the sign changes of DeltaPsi ----------------------
dps = @(l) dPsi(l, cf, m, gk, t);
NS   = 12*max(m,4);
lsc  = lam_of_th(linspace(1e-6, pi-1e-6, NS));
vsc  = arrayfun(dps, lsc);
sgn  = find(vsc(1:end-1).*vsc(2:end) < 0);
cross = zeros(1, numel(sgn));
for i = 1:numel(sgn)
    % Brent (fzero) on the bracketing pair: ~10 evaluations instead of 60
    cross(i) = fzero(dps, [lsc(sgn(i)), lsc(sgn(i)+1)]);
end

% ---------- 4. Chebyshev-mapped, panel-split Gauss quadrature -----------
[xg, wg] = gauss_legendre(ng);
edges = unique([0, th_of_lam(cross), pi]);
c_inf = 0;  nev = numel(vsc) + 12*numel(cross);
for i = 1:numel(edges)-1
    sub = linspace(edges(i), edges(i+1), npan+1);
    for j = 1:npan
        c = (sub(j)+sub(j+1))/2;  h = (sub(j+1)-sub(j))/2;
        tt = c + h*xg;
        vv = arrayfun(@(s) abs(dps(lam_of_th(s))), tt);
        c_inf = c_inf + h*sum(wg.*vv.*sin(tt));
        nev = nev + ng;
    end
end
c_inf = c_inf*hw/pi;

% ---------- 5. internal error estimate: bandwidth sensitivity -----------
t2   = max(8, t-6);
dps2 = @(l) dPsi(l, cf, m, gk, t2);
ls   = lam_of_th(linspace(0.05, pi-0.05, 40));
tstab = max(abs(arrayfun(dps,ls) - arrayfun(dps2,ls)));
if verbose
    fprintf('crossings: %d,  integrand evaluations: %d\n', numel(cross), nev);
    fprintf('bandwidth sensitivity (t=%d vs %d): %.1e\n', t, t2, tstab);
    fprintf('c_inf = %.8e     bracket [%.3e, %.3e]\n', c_inf, floorb, ...
            (B-A)/pi*max(abs(vsc)));
end

% ---------- 6. optional cross-check at moderate n ----------------------
ds = NaN;
if docheck
    n  = 40;
    Tf = toeplitz([a0, apos, zeros(1,n-m-1)], [a0, aneg, zeros(1,n-m-1)]);
    kg = min(n-1, K);
    Tg = toeplitz([gk(1), gk(2:kg+1), zeros(1,n-kg-1)]);
    ds = sum(abs(sort(real(eig(Tf))) - sort(eig(Tg))));
    if verbose, fprintf('check: d_sigma(n=40) = %.4e\n', ds); end
end

out = struct('c_inf',c_inf,'Delta2',Delta2,'floor',floorb,'A',A,'B',B, ...
             'gk',gk,'t',t,'crossings',cross,'nevals',nev, ...
             'tstab',tstab,'ds_check',ds);
end

% ======================================================================
function d = dPsi(lam, cf, m, gk, t)
% Phase mismatch Psi_g - Psi_f, in cancellation-free form (the leading
% alpha of each phase cancels, so it never enters).

% --- f side: roots of z^m (f_m - lam), degree 2m -----------------------
c = cf;  c(m+1) = c(m+1) - lam;
r = roots(c);  [~, ix] = sort(abs(r));  r = r(ix);
if imag(r(m)) > 0, zc = r(m); else, zc = r(m+1); end
zb = conj(zc);
inner = r(1:m-1);  outer = r(m+2:end);
Sf = 0.5*sum(angle( zc*(zb - inner) ./ (zb*(zc - inner)) )) ...
   - 0.5*sum(angle( (outer - zb) ./ (outer - zc) ));

% --- g side: palindromic => degree-t colleague matrix in Chebyshev basis
a = [gk(1) - lam, 2*gk(2:t+1)];              % Chebyshev coeffs, u = 2x
C = zeros(t);
C(1,2) = 1;
for k = 2:t-1, C(k,k-1) = 0.5;  C(k,k+1) = 0.5;  end
C(t,t-1) = 0.5;
C(t,:) = C(t,:) - a(1:t)/(2*a(t+1));
u = 2*eig(C);
s = sqrt(u.^2 - 4);
z = (u - s)/2;  zo = (u + s)/2;
sw = abs(z) > 1;  z(sw) = zo(sw);            % keep the root inside the disc
[~, j] = min(abs(abs(z) - 1));               % the unimodular central root
zg = z(j);  if imag(zg) < 0, zg = conj(zg); end
z(j) = [];
Sg = sum(angle( zg*(1 - zg*z) ./ (zg - z) ));

d = Sg - Sf;
end

function z = central_root(lam, cf, m)
c = cf;  c(m+1) = c(m+1) - lam;
r = roots(c);  [~, ix] = sort(abs(r));  r = r(ix);
if imag(r(m)) > 0, z = r(m); else, z = r(m+1); end
end

function [x, w] = gauss_legendre(n)
k = 1:n-1;  b = k./sqrt(4*k.^2 - 1);
[V, D] = eig(diag(b,1) + diag(b,-1));
[x, i] = sort(diag(D));  x = x(:).';  w = 2*V(1,i).^2;
end

function y = edge_of(y0, A, B)
if abs(y0-A) < abs(y0-B), y = A; else, y = B; end
end

function [A, B] = locate_band(cf, m)
scan = linspace(0.3, 2.0, 4000);
inb  = false(size(scan));
for i = 1:numel(scan), inb(i) = in_band(scan(i), cf, m); end
ix = find(inb);
assert(~isempty(ix), 'no band found: widen the scan interval');
A = bisect(scan(ix(1)),   scan(max(ix(1)-1,1)), cf, m);
B = bisect(scan(ix(end)), scan(min(ix(end)+1,numel(scan))), cf, m);
end

function tf = in_band(lam, cf, m)
c = cf;  c(m+1) = c(m+1) - lam;
r = roots(c);  [~, ix] = sort(abs(r));  r = r(ix);
tf = abs(imag(r(m))) > 1e-12*abs(r(m)) && ...
     abs(abs(r(m)) - abs(r(m+1))) < 1e-9*abs(r(m));
end

function x = bisect(xin, xout, cf, m)
for it = 1:80
    xm = (xin + xout)/2;
    if in_band(xm, cf, m), xin = xm; else, xout = xm; end
end
x = xin;
end
