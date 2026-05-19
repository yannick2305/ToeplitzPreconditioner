function r = sampleUniformlyVec(x, y, m, k)

d  = 2*k + 1;
w  = d + 1;          % stencil width = 2k+2
n  = length(x);

% --- Polar conversion (vectorized) ---
rho   = hypot(x, y);
theta = acos(x ./ rho);

% --- Extend data at both ends ---
rho2   = [rho(k+1:-1:1),    rho,    rho(n:-1:n-k)  ];
theta2 = [-theta(k+1:-1:1), theta,  2*pi-theta(n:-1:n-k)];
n2     = length(theta2);     % = n + w

% --- Query angles ---
u = (0:m-1) * pi / (m-1);   % 1 x m

% --- Locate stencil for every u simultaneously ---
% histcounts gives mu(i) s.t. theta2(mu(i)) <= u(i) < theta2(mu(i)+1)
[~, ~, mu] = histcounts(u, theta2);
mu = min(mu, n2 - k - 1);   % clamp: last point u=pi may hit right edge

% --- Gather stencil nodes and values: m x w matrices ---
idx   = mu(:) + (-k : k+1); % m x w  (each row is one stencil)
nodes = theta2(idx);          % m x w
vals  = rho2(idx);            % m x w

% --- Vectorized Lagrange interpolation ---
u_col = u(:);                 % m x 1

% Differences (u_i - node_{i,l}) for all i, l:  m x w
du = u_col - nodes;

% Numerator of basis l: full product then divide out l-th factor
full_prod = prod(du, 2);      % m x 1
num       = full_prod ./ du;  % m x w

% Denominator: prod_{l' != l} (node_{i,l} - node_{i,l'})
% Build m x w x w difference tensor using implicit broadcasting
D = reshape(nodes, [m,w,1]) - reshape(nodes, [m,1,w]);  % m x w x w
D = D + reshape(eye(w, 'logical'), [1,w,w]);             % diagonal 0 -> 1
denom = prod(D, 3);           % m x w

r = sum((num ./ denom) .* vals, 2)';   % 1 x m