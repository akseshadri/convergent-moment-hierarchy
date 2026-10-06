function mu = truncvar(mi, Sigma, alpha, n1, n2)
% TRUNCVAR  Variance approximation mu^(n1,n2) for uniform weights.
%
%   mu = truncvar(mi, Sigma, alpha, n1, n2) returns the truncated variance
%   approximation mu_2^(n1,n2) = E[T_{n1,n2}^2] for uniform weights under MCAR
%   reporting with availability alpha, for any truncation orders n1, n2 >= 0.
%   n1 counts the mixed (numerator-denominator) terms and n2 the denominator
%   terms retained in the centered ratio expansion.
%
%   With uniform weights the denominator fluctuation V = K/N - alpha depends
%   only on the number K of reporting sites, and conditional on K = k the
%   numerator moments follow from sampling k of N sites without replacement.
%   The truncation is therefore evaluated exactly as a sum over k = 0..N.
%   Requires N >= 2.

if alpha == 1                             % complete reporting: all truncations equal Var(R)
    mu = sum(Sigma(:))/numel(mi)^2; return
end
N   = numel(mi); mi = mi(:);
Sd  = trace(Sigma); So = sum(Sigma(:)) - Sd;
sm  = sum(mi); sm2 = sum(mi.^2); mr = sm/N;
k   = (0:N)';
pk  = exp(gammaln(N+1) - gammaln(k+1) - gammaln(N-k+1) + k*log(alpha) + (N-k)*log(1-alpha));
mS  = alpha; mR = alpha*mr;               % means of S and R
V   = k/N - alpha;                        % denominator fluctuation given K = k
ER  = k*mr/N;                             % E[R | K = k]
ER2 = ((k/N)*(Sd + sm2) + (k.*(k-1)/(N*(N-1)))*(So + sm^2 - sm2))/N^2;   % E[R^2 | K = k]
EU  = ER - mR;                            % E[U | K = k], U = R - mR
EU2 = ER2 - 2*mR*ER + mR^2;               % E[U^2 | K = k]

% truncated centered series: T = c(K)*U + d(K)
c = zeros(N+1,1); d = zeros(N+1,1);
for n = 0:n1                              % mixed terms (-1)^n (U V^n - E[U V^n]) / mS^(n+1)
    c = c + (-1)^n * V.^n / mS^(n+1);
    d = d - (-1)^n * sum(pk.*EU.*V.^n) / mS^(n+1);
end
for n = 1:n2                              % denominator terms (-1)^n mR (V^n - E[V^n]) / mS^(n+1)
    d = d + (-1)^n * mR * (V.^n - sum(pk.*V.^n)) / mS^(n+1);
end
mu = sum(pk .* (c.^2.*EU2 + 2*c.*d.*EU + d.^2));
