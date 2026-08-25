function T = truncations_general(beta, mi, sigmamat, alpha)
% TRUNCATIONS_GENERAL  Arbitrary-weight MAR closed forms for the weighted
%   ratio estimator: the (0,1),(1,1),(0,2) variance truncations and rho.
%
%   T = TRUNCATIONS_GENERAL(beta, mi, sigmamat, alpha) evaluates the closed
%   forms for general weights beta and an arbitrary field covariance sigmamat;
%   it reduces to Scenario I for uniform weights. kappa = Var(UV) is assembled
%   from the coefficient blocks T_ij = E[s_i s_j S] and Q_ij = E[s_i s_j S^2].
%
%   Inputs:
%     beta      — [N x 1] weights (uniform is ones(N,1)/N).
%     mi        — [N x 1] site means m_i = E[r_i].
%     sigmamat  — [N x N] field covariance; diag sigma_i^2, off-diag sigma_ij.
%     alpha     — reporting (availability) probability, in (0,1).
%
%   Outputs (fields of struct T):
%     mu00,mu01,mu11,mu02  — the (0,0),(0,1),(1,1),(0,2) variance truncations.
%     rho                  — relative denominator fluctuation sqrt(Var(S))/E[S].
%
%   Verified against exact 2^N enumeration.
    beta = beta(:); mi = mi(:); N = numel(beta);
    sigmai2 = diag(sigmamat);
    m2i = sigmai2 + mi.^2;                         % E[r_i^2]
    m2  = sigmamat + mi*mi.'; m2(1:N+1:end) = m2i; % E[r_i r_j]

    W2 = sum(beta.^2); W3 = sum(beta.^3); W4 = sum(beta.^4);
    m_r = sum(beta.*mi); M2 = sum(beta.^2.*mi); M3 = sum(beta.^3.*mi);
    H3 = sum(beta.^3.*mi.^2); V3 = sum(beta.^3.*sigmai2);

    B = beta*beta.'; Boff = B; Boff(1:N+1:end) = 0;
    Csig = sum(sum(Boff .* (beta+beta.') .* sigmamat));

    m_S = alpha; m_R = alpha*m_r;
    Var_S = alpha*(1-alpha)*W2;
    Cov_RS = alpha*(1-alpha)*M2;
    Var_R = alpha*sum(sigmai2.*beta.^2) + alpha^2*sum(sum(Boff.*sigmamat)) + alpha*(1-alpha)*sum(beta.^2.*mi.^2);

    T.mu00 = Var_R/m_S^2;
    T.mu01 = T.mu00 + m_R^2/m_S^4*Var_S - 2*m_R/m_S^3*Cov_RS;

    mu_1_2 = alpha*(1-alpha)*(1-2*alpha)*M3;
    mu_2_1 = alpha*(1-alpha)*(V3 + alpha*Csig + (1-2*alpha)*H3);

    bi = beta*ones(1,N); bj = ones(N,1)*beta.';
    Tblk = alpha^3 + alpha^2*(1-alpha)*(bi+bj);
    Qblk = alpha^2*(bi+bj).^2 + alpha^3*(2*bi+2*bj+W2-3*bi.^2-3*bj.^2-4*bi.*bj) ...
         + alpha^4*(1-2*bi-2*bj-W2+2*bi.^2+2*bj.^2+2*bi.*bj);
    Tdiag = alpha^2 + alpha*(1-alpha)*beta;
    Qdiag = alpha*beta.^2 + alpha^2*(2*beta+W2-3*beta.^2) + alpha^3*(1-2*beta+2*beta.^2-W2);
    Tblk(1:N+1:end) = Tdiag; Qblk(1:N+1:end) = Qdiag;

    bm = beta.*mi;
    E_RS   = sum(bm.*diag(Tblk));   E_RS2  = sum(bm.*diag(Qblk));
    E_R2S  = sum(sum(B.*m2.*Tblk)); E_R2S2 = sum(sum(B.*m2.*Qblk));
    E_S2   = alpha^2 + alpha*(1-alpha)*W2;
    E_R2   = alpha*sum(beta.^2.*m2i) + alpha^2*sum(sum(Boff.*m2));
    E_U2V2 = E_R2S2 - 2*m_R*E_RS2 + m_R^2*E_S2 - 2*m_S*E_R2S + 4*m_R*m_S*E_RS + m_S^2*E_R2 - 3*m_R^2*m_S^2;
    kappa = E_U2V2 - Cov_RS^2;

    T.mu11 = T.mu01 - 2/m_S^3*mu_2_1 + 2*m_R/m_S^4*mu_1_2 + 1/m_S^4*kappa;

    E_V3 = alpha*(1-alpha)*(1-2*alpha)*W3;
    E_V4 = alpha*(1-alpha)*(1-3*alpha*(1-alpha))*W4 + 3*alpha^2*(1-alpha)^2*(W2^2-W4);
    Var_V2 = E_V4 - (alpha*(1-alpha)*W2)^2;
    T.mu02 = T.mu01 + 2*m_R/m_S^4*mu_1_2 - 2*m_R^2/m_S^5*E_V3 + m_R^2/m_S^6*Var_V2;

    T.rho = sqrt((1-alpha)/alpha*W2);
end
