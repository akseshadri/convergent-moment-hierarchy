function v = varexact(mi, Sigma, alpha)
% VAREXACT  Exact variance of the spatial average for uniform weights.
%
%   v = varexact(mi, Sigma, alpha) returns Var(rhat) for uniform weights
%   beta_i = 1/N under MCAR reporting with availability alpha, using the
%   convention rhat = mean(mi) when no site reports.  mi is the N x 1 vector
%   of site means and Sigma the N x N covariance matrix.
%
%   Conditional on K = k reporting sites, rhat is the mean of k sites drawn
%   without replacement, so the variance is a single sum over k (main text,
%   Eq. for the exact variance with uniform weights).

if alpha == 1                             % complete reporting
    v = sum(Sigma(:))/numel(mi)^2; return
end
N  = numel(mi);
Sd = trace(Sigma);                        % sum of site variances
So = sum(Sigma(:)) - Sd;                  % sum of off-diagonal covariances
Vm = mean((mi - mean(mi)).^2);            % spatial variance of site means
k  = (1:N)';
pk = exp(gammaln(N+1) - gammaln(k+1) - gammaln(N-k+1) + k*log(alpha) + (N-k)*log(1-alpha));
v  = sum(pk .* (Sd./(k*N) + (k-1)./(k*N*(N-1))*So + (N-k)./((N-1)*k)*Vm));
