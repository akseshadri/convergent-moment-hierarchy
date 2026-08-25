function T = truncations_scenarioI(alpha, N, Sum_sigma_diag, Sum_sigma_offdiag, M2_m, Var_m, m_r)
% TRUNCATIONS_SCENARIOI  Uniform-weight (Scenario I) truncated variance
%   approximations for the weighted ratio estimator under MAR.
%
%   T = TRUNCATIONS_SCENARIOI(alpha, N, Sum_sigma_diag, Sum_sigma_offdiag, ...
%                             M2_m, Var_m, m_r) evaluates the closed-form
%   variance truncations from the field-moment sums, for uniform weights
%   beta_i = 1/N and reporting probability alpha. The homogeneous ("Scenario
%   III") case is the special input Sum_sigma_diag = N*sigma^2,
%   Sum_sigma_offdiag = 0, M2_m = N*m^2, Var_m = 0, m_r = m.
%
%   Inputs:
%     alpha              — reporting (availability) probability, in (0,1).
%     N                  — number of sites.
%     Sum_sigma_diag     — sum_i sigma_i^2 (sum of site variances).
%     Sum_sigma_offdiag  — sum_(i~=j) sigma_ij (off-diagonal covariances).
%     M2_m               — sum_i m_i^2 (sum of squared site means).
%     Var_m              — mean_i (m_i - m_r)^2 (spatial variance of means).
%     m_r                — mean_i m_i (grand mean).
%
%   Outputs (fields of struct T):
%     mu00,mu01,mu10,mu11  — the (0,0),(0,1),(1,0),(1,1) variance truncations.
%     mu12,mu21,kappa      — mixed moments mu_(1,2)=E[U V^2], mu_(2,1)=E[U^2 V],
%                          kappa=Var(UV), with U=R-E[R], V=S-E[S].
%     t00,t01,t10,t11      — the constituent terms of each truncation (the
%                          vectors plotted by the decomposition figures).
%
%   Verified against exact 2^N enumeration.
    m_R = alpha*m_r;

    % (0,0): baseline variance + spatial covariance + average squared mean
    t1  = 1/(alpha*N^2)*Sum_sigma_diag;
    t2  = (1/N^2)*Sum_sigma_offdiag;
    t3  = ((1-alpha)/alpha)*(1/N^2)*M2_m;
    T.mu00 = t1 + t2 + t3;
    T.t00  = [t1 t2 t3];

    % (0,1): replaces the average squared mean by the mean heterogeneity
    t3b = ((1-alpha)/alpha)*(1/N)*Var_m;
    T.mu01 = t1 + t2 + t3b;
    T.t01  = [t1 t2 t3b];

    % mixed moments (eqs. 3.11, 3.12, 3.15)
    T.mu12  = (alpha*(1-alpha)*(1-2*alpha)/N^2)*m_r;
    T.mu21  = (alpha*(1-alpha)/N^3)*(Sum_sigma_diag + 2*alpha*Sum_sigma_offdiag + (1-2*alpha)*M2_m);
    T.kappa = (alpha*(1-alpha)/N^4)*((1+alpha*(N-2))*Sum_sigma_diag ...
              + (4*alpha+alpha^2*(N-6))*Sum_sigma_offdiag ...
              + (1+alpha*(1-alpha)*(N-6))*M2_m + alpha*(1-alpha)*N^2*m_r^2);

    % (1,0): (0,0) plus the leading coupling corrections
    c4 = -2/alpha^3*T.mu21;
    c5 =  1/alpha^4*T.kappa;
    T.mu10 = T.mu00 + c4 + c5;
    T.t10  = [t1 t2 t3 c4 c5];

    % (1,1): (0,1) plus the coupling corrections
    c5b = 2/alpha^4*m_R*T.mu12;
    T.mu11 = T.mu01 + c4 + c5b + c5;
    T.t11  = [t1 t2 t3b c4 c5b c5];
end
