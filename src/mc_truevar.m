function [varmc, mu21mc, mu12mc, kappamc] = mc_truevar(datause, beta, alpha, Kensemb)
% MC_TRUEVAR  Monte-Carlo moments of rhat = R/S under MAR sampling.
%
%   [varmc, mu21mc, mu12mc, kappamc] = MC_TRUEVAR(datause, beta, alpha, Kensemb)
%   holds the field DATAUSE fixed and resamples only the reporting pattern
%   s_ij ~ Bernoulli(alpha). On an all-missing day (S = 0) the convention
%   rhat = m_r is applied, keeping rhat bounded; that event is exponentially
%   rare. Averages are taken over Kensemb ensemble members.
%
%   Inputs:
%     datause  — [Tuse x N] field (days x sites); held fixed.
%     beta     — [N x 1] weights (uniform or general).
%     alpha    — reporting (availability) probability, in (0,1).
%     Kensemb  — number of Monte-Carlo ensemble members.
%
%   Outputs:
%     varmc                    — expected temporal variance of rhat.
%     mu21mc, mu12mc, kappamc  — MC mixed moments mu_(2,1)=E[U^2 V],
%                              mu_(1,2)=E[U V^2], kappa=Var(UV), with
%                              U=R-E[R], V=S-E[S] (optional outputs).
    [Tuse, N] = size(datause);
    beta = beta(:);
    betamat = repmat(beta', [Tuse 1]);
    m_r = sum(beta .* mean(datause,1)');
    varmclist = NaN(Kensemb,1);
    mu21list  = NaN(Kensemb,1); mu12list = NaN(Kensemb,1); kappalist = NaN(Kensemb,1);
    parfor i = 1:Kensemb
        sij   = binornd(1, alpha, Tuse, N);
        Rjsum = sum(betamat .* datause .* sij, 2);
        Sjsum = sum(betamat .* sij, 2);
        rlist = Rjsum ./ Sjsum;
        rlist(Sjsum==0) = m_r;
        varmclist(i) = var(rlist);
        Ulist = Rjsum - mean(Rjsum);
        Vlist = Sjsum - mean(Sjsum);
        mu21list(i) = mean(Ulist.^2 .* Vlist);
        mu12list(i) = mean(Ulist .* Vlist.^2);
        kappalist(i) = mean(Ulist.^2 .* Vlist.^2) - (mean(Ulist.*Vlist))^2;
    end
    varmc  = mean(varmclist);
    mu21mc = mean(mu21list);
    mu12mc = mean(mu12list);
    kappamc = mean(kappalist);
end
