function v = ratio_var_exact_uniform(m, sig2, alpha)
% RATIO_VAR_EXACT_UNIFORM  Exact Var(rhat) for uniform weights and an
%   independent field, with no Monte Carlo.
%
%   v = RATIO_VAR_EXACT_UNIFORM(m, sig2, alpha) returns the exact variance of
%   rhat = R/S with beta_i = 1/N and the convention rhat = mean(m) when S = 0.
%   It conditions on the number k of reporting sites: for uniform weights
%   S = k/N is fixed given |A| = k, and averaging over the C(N,k) equally
%   likely reporting sets gives, for k >= 1,
%     E[rhat | k]   = mbar,
%     E[rhat^2 | k] = mbar^2 + (N-k)/((N-1)*k)*sm2 + sb2/k,
%   with sm2 the spatial variance of the site means and sb2 their mean
%   variance. Summing against the Binomial(N,alpha) law of k (the k = 0 term
%   contributes 0) yields Var(rhat).
%
%   Inputs:
%     m      — [N x 1] site means m_i = E[r_i].
%     sig2   — [N x 1] site variances sigma_i^2 (independent field).
%     alpha  — reporting (availability) probability, in (0,1).
%
%   Output:
%     v      — exact Var(rhat).
%
%   Verified against brute-force 2^N enumeration.
    m = m(:); sig2 = sig2(:); N = numel(m);
    mbar = mean(m);
    sm2  = mean((m - mbar).^2);          % spatial variance of the site means
    sb2  = mean(sig2);                   % mean site variance
    k    = (1:N).';
    % Binomial pmf via log-gamma (stable for large N, no toolbox dependence)
    logpk = gammaln(N+1) - gammaln(k+1) - gammaln(N-k+1) ...
            + k*log(alpha) + (N-k)*log(1-alpha);
    pk   = exp(logpk);
    term = (N - k)./((N - 1).*k)*sm2 + sb2./k;
    v    = sum(pk .* term);              % k = 0 contributes 0 (rhat = mbar there)
end
