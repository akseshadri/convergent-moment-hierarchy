function [v, mu, mm] = enumvar(beta, mi, Sigma, alpha, orders)
% ENUMVAR  Exact variance and truncations by enumerating reporting patterns.
%
%   [v, mu, mm] = enumvar(beta, mi, Sigma, alpha, orders) averages over all
%   2^N reporting patterns for arbitrary weights beta (small N only).
%     v      exact Var(rhat), with rhat = sum(beta.*mi) when no site reports
%     mu     truncations mu^(n1,n2) for each row [n1 n2] of orders
%     mm     exact mixed moments [mu12 mu21 kappa]

N = numel(beta); beta = beta(:); mi = mi(:);
P = dec2bin(0:2^N-1) - '0';               % all reporting patterns (2^N x N)
kk = sum(P,2);
p  = alpha.^kk .* (1-alpha).^(N-kk);      % pattern probabilities
Q  = P .* repmat(beta',size(P,1),1);      % beta_i s_i
S  = sum(Q,2);
mS = alpha; mR = alpha*(beta'*mi);
V  = S - mS;
EU  = Q*mi - mR;                          % E[U | pattern]
EU2 = sum((Q*Sigma).*Q,2) + EU.^2;        % E[U^2 | pattern]

% exact variance of rhat
cm = repmat(beta'*mi, size(S)); cv = zeros(size(S));
on = S > 0;
cm(on) = (Q(on,:)*mi)./S(on);
cv(on) = sum((Q(on,:)*Sigma).*Q(on,:),2)./S(on).^2;
v = sum(p.*(cv + cm.^2)) - sum(p.*cm)^2;

% truncations
mu = zeros(size(orders,1),1);
for r = 1:size(orders,1)
    c = zeros(size(V)); d = zeros(size(V));
    for n = 0:orders(r,1)
        c = c + (-1)^n * V.^n / mS^(n+1);
        d = d - (-1)^n * sum(p.*EU.*V.^n) / mS^(n+1);
    end
    for n = 1:orders(r,2)
        d = d + (-1)^n * mR * (V.^n - sum(p.*V.^n)) / mS^(n+1);
    end
    mu(r) = sum(p.*(c.^2.*EU2 + 2*c.*d.*EU + d.^2));
end

% mixed moments
EUV = sum(p.*EU.*V);
mm = [sum(p.*EU.*V.^2), sum(p.*EU2.*V), sum(p.*EU2.*V.^2) - EUV^2];
