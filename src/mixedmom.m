function mm = mixedmom(beta, mi, Sigma, alpha)
% MIXEDMOM  Closed-form MCAR mixed moments [mu12 mu21 kappa].
%
%   mm = mixedmom(beta, mi, Sigma, alpha) evaluates the closed forms of
%   mu12 = E[U V^2], mu21 = E[U^2 V], and kappa = Var(U V) for arbitrary
%   weights under MCAR (main text, Section 3; Supplementary Information S1).

beta = beta(:); mi = mi(:); N = numel(beta);
W2 = sum(beta.^2); mr = beta'*mi;
M2 = sum(beta.^2.*mi); M3 = sum(beta.^3.*mi);
H3 = sum(beta.^3.*mi.^2); V3 = sum(beta.^3.*diag(Sigma));
off = ~eye(N);
Bsum = repmat(beta,1,N) + repmat(beta',N,1);
Csig = sum(sum(off.*(beta*beta').*Bsum.*Sigma));
mu12 = alpha*(1-alpha)*(1-2*alpha)*M3;
mu21 = alpha*(1-alpha)*(V3 + alpha*Csig + (1-2*alpha)*H3);

% kappa from the reporting-coefficient blocks T_ij and Q_ij
m2 = Sigma + mi*mi';
bi = repmat(beta,1,N); bj = repmat(beta',N,1);
Tblk = alpha^3 + alpha^2*(1-alpha)*(bi + bj);
Tblk(1:N+1:end) = alpha^2 + alpha*(1-alpha)*beta;
Qblk = alpha^2*(bi+bj).^2 + alpha^3*(2*(bi+bj).*(1-bi-bj) + W2 - bi.^2 - bj.^2) ...
     + alpha^4*((1-bi-bj).^2 - (W2 - bi.^2 - bj.^2));
Qblk(1:N+1:end) = alpha*beta.^2 + alpha^2*(2*beta + W2 - 3*beta.^2) ...
     + alpha^3*(1 - 2*beta + 2*beta.^2 - W2);
mS = alpha; mR = alpha*mr;
ES2 = alpha*W2 + alpha^2*(1 - W2);
ER2 = alpha*sum(beta.^2.*diag(m2)) + alpha^2*(beta'*m2*beta - sum(beta.^2.*diag(m2)));
ERS = sum(beta.*mi.*diag(Tblk)); ERS2 = sum(beta.*mi.*diag(Qblk));
ER2S = sum(sum((beta*beta').*m2.*Tblk)); ER2S2 = sum(sum((beta*beta').*m2.*Qblk));
EU2V2 = ER2S2 - 2*mS*ER2S + mS^2*ER2 - 2*mR*ERS2 + 4*mR*mS*ERS + mR^2*ES2 - 3*mR^2*mS^2;
covRS = alpha*(1-alpha)*M2;
mm = [mu12, mu21, EU2V2 - covRS^2];
