%% PFigs_Section2.m — Section 2 figures: rho governs convergence and accuracy
%
% Companion code for:
%   Seshadri, A. K. and Pal Majumder, A. (2026). A convergent moment
%   hierarchy for the variance of weighted spatial-mean estimators under
%   missing-at-random sampling.
%
% Produces the two figures illustrating Section-2 results showing that 
% a single dimensionless parameter, the relative denominator fluctuation
%   rho = sqrt(Var(S))/E[S] = sqrt((1-alpha)/(alpha*N))   (uniform weights),
% governs both the convergence of the moment expansion and the accuracy of
% its low-order truncations.
%
%   Figure 1 (rho master):
%     (a) rho-collapse of the (0,0) truncation error onto m^2*rho^2 for all
%         (alpha, N);
%     (b) the sharp rates: (0,0) error ~ O(rho^2), (0,1) error ~ O(rho^4);
%     (c) convergence ~ 1/N at every fixed alpha.
%   Figure 2 (good / rare split):
%     (a) the denominator fluctuation (S-ES)/ES with the good event
%         |(S-ES)/ES| <= tau shaded;
%     (b) the exact rare-event probability P(G^c) below the Chebyshev
%         (rho^2/tau^2) and Bernstein bounds.
%
% Section 2's claims are asymptotic-in-N scaling results, clearest in the
% homogeneous ("Scenario III") setting. We drive them with a homogeneous field
% whose scale sigma is taken from the IMD rainfall data, and evaluate same
% closed forms as the other scripts through TRUNCATIONS_SCENARIOI with the
% homogeneous moment sums (Sum_sigma_diag = N*sigma^2, Sum_sigma_offdiag = 0,
% M2_m = N*m^2, Var_m = 0, m_r = m). The (0,1) error is O(rho^4), below the
% Monte-Carlo floor, so the rho^2 vs rho^4 rates are shown with analytic
% increments; Monte Carlo confirms the resolvable (0,0) error. P(G^c) is the
% exact binomial tail.
%
% Uses: truncations_scenarioI, mc_truevar, apply_fig_style, panel_label.
%
% Data: 1 x 1 deg IMD gridded daily rainfall (../data/indiadat.mat), used only
% to set the rainfall scale sigma.
%
% Requirements: Statistics and Machine Learning Toolbox (binornd, binopdf);
% Parallel Computing Toolbox optional.


clear                                   

lw = 1.5; ms = 8;                       % consistent line width and marker size
cMC = 'k'; c00 = 'r'; c01 = 'b'; c11 = 'c';   % consistent colours

%% Load data (used only to set the rainfall scale sigma)
inddat = load('../data/indiadat.mat');
datmat = NaN(111,357,365);
for i = 1901:2011
    datmat(i-1900,:,:) = single(cell2mat(inddat.indiarainmodel(i)));
end
data = reshape(datmat,[357 111*365])';   % [time,location]
% (portable alternative:  data = readmatrix('data.txt');)

Tuse = 365;                              % number of daily fields used (one year)
sigmai2 = var(data(1:Tuse,:),[],1)*(Tuse-1)/Tuse;
sigma = sqrt(mean(sigmai2));
m = sigma;                               % "mean comparable to standard deviation"
fprintf('Rainfall-derived sigma = %.3f; homogeneous field uses m = sigma = %.3f\n', sigma, m);

%% ================= FIGURE 1 : rho governs =================
figure(1); clf

% ---- (a) rho-collapse of the (0,0) error across (alpha,N) ----
subplot(1,3,1)
Nlist = [40 80 160 320];
alphalist = [0.15 0.25 0.35 0.5 0.65 0.8 0.9];
markers = {'o','s','d','^'};
for nn = 1:numel(Nlist)
    N = Nlist(nn);
    rho_list = NaN(numel(alphalist),1); err_list = NaN(numel(alphalist),1);
    for k = 1:numel(alphalist)
        alpha = alphalist(k);
        T = truncations_scenarioI(alpha, N, N*sigma^2, 0, N*m^2, 0, m);
        rho_list(k) = sqrt((1-alpha)/(alpha*N));
        err_list(k) = abs(T.mu00 - T.mu01);
    end
    loglog(rho_list, err_list, markers{nn}, 'MarkerSize', ms, 'LineWidth', 1), hold on
end
rr = logspace(-1.7,-0.3,50);
loglog(rr, m^2*rr.^2, 'k--', 'LineWidth', lw)
xlabel('\rho = sqrt((1-\alpha)/(\alpha N))'), ylabel('|Var - \mu_2^{(0,0)}|')
legend('N=40','N=80','N=160','N=320','m^2\rho^2')
panel_label('(a)'), apply_fig_style(gca)

% ---- (b) sharp rates: rho via N at fixed alpha = 0.3 ----
subplot(1,3,2)
alpha = 0.3; Nlist = [20 30 45 70 100 150 220 300 357]; Kensemb = 6000;
rho_list = NaN(numel(Nlist),1); incr01 = rho_list; incr11 = rho_list; errMC = rho_list;
for nn = 1:numel(Nlist)
    N = Nlist(nn);
    T = truncations_scenarioI(alpha, N, N*sigma^2, 0, N*m^2, 0, m);
    rho_list(nn) = sqrt((1-alpha)/(alpha*N));
    incr01(nn) = abs(T.mu00 - T.mu01);
    incr11(nn) = abs(T.mu11 - T.mu01);
    datause = m + sigma*randn(Tuse,N);
    errMC(nn) = abs(mc_truevar(datause, ones(N,1)/N, alpha, Kensemb) - T.mu00);
end
loglog(rho_list, incr01, [c00 '-o'], 'LineWidth', lw, 'MarkerSize', 6), hold on
loglog(rho_list, errMC, [cMC '+'], 'LineWidth', 2, 'MarkerSize', 11)
loglog(rho_list, incr11, [c01 '-s'], 'LineWidth', lw, 'MarkerSize', 6)
loglog(rho_list, m^2*rho_list.^2, [c00 '--'])
loglog(rho_list, incr11(end)*(rho_list/rho_list(end)).^4, [c01 '--'])
xlabel('\rho (via N, \alpha = 0.3)'), ylabel('increment / error')
legend('|\mu_2^{(0,0)}-\mu_2^{(0,1)}|','|Var_{MC}-\mu_2^{(0,0)}|','|\mu_2^{(1,1)}-\mu_2^{(0,1)}|','slope 2','slope 4')
panel_label('(b)'), apply_fig_style(gca)

% ---- (c) convergence: error ~ 1/N at every fixed alpha ----
subplot(1,3,3)
Nlist = [15 25 40 60 90 140 220 320]; alphalist = [0.1 0.3 0.5 0.9];
alcol = {c00, c01, [0 0.6 0], c11}; Kensemb = 5000;
for k = 1:numel(alphalist)
    alpha = alphalist(k);
    ana = NaN(numel(Nlist),1); mc = ana;
    for nn = 1:numel(Nlist)
        N = Nlist(nn);
        T = truncations_scenarioI(alpha, N, N*sigma^2, 0, N*m^2, 0, m);
        ana(nn) = abs(T.mu00 - T.mu01);
        datause = m + sigma*randn(Tuse,N);
        mc(nn) = abs(mc_truevar(datause, ones(N,1)/N, alpha, Kensemb) - T.mu00);
    end
    loglog(Nlist, ana, '-', 'Color', alcol{k}, 'LineWidth', lw), hold on
    loglog(Nlist, mc, 'o', 'Color', alcol{k}, 'MarkerSize', 5)
end
xlabel('N'), ylabel('|Var - \mu_2^{(0,0)}|')
legend('\alpha=0.1','','\alpha=0.3','','\alpha=0.5','','\alpha=0.9','')
panel_label('(c)'), apply_fig_style(gca)


%% ================= FIGURE 2 : good / rare split =================
%% Exact binomial tail, N*S ~ Binomial(N,alpha), good event {|(S-ES)/ES|<=tau}.
tau = 0.5;
figure(2); clf

% ---- (a) distribution of the denominator fluctuation ----
subplot(1,2,1)
Nlist = [50 200 357]; ncol = {c00, c01, [0 0.6 0]};
for nn = 1:numel(Nlist)
    N = Nlist(nn); alpha = 0.5;
    k = 0:N; pmf = binopdf(k,N,alpha);
    x = (k/N - alpha)/alpha;
    dens = pmf/((1/N)/alpha);
    plot(x, dens, '-', 'Color', ncol{nn}, 'LineWidth', lw), hold on
end
yl = ylim; patch([-tau -tau tau tau], [0 yl(2) yl(2) 0], 'k', 'FaceAlpha', 0.06, 'EdgeColor', 'none')
xline(tau, ':k'), xline(-tau, ':k')
xlim([-1.2 1.2]), xlabel('(S - ES) / ES'), ylabel('density (exact)')
legend('N=50','N=200','N=357')
panel_label('(a)'), apply_fig_style(gca)

% ---- (b) exact P(G^c) vs the Chebyshev and Bernstein bounds ----
subplot(1,2,2)
Nlist = [25 50 100 200 357 600 1000]; alpha = 0.5;
PGc = NaN(numel(Nlist),1); Cheb = PGc; Bern = PGc;
for nn = 1:numel(Nlist)
    N = Nlist(nn);
    k = 0:N; pmf = binopdf(k,N,alpha); Sfrac = k/N;
    PGc(nn) = sum(pmf(abs(Sfrac-alpha) > tau*alpha));
    rho2 = (1-alpha)/(alpha*N);
    Cheb(nn) = rho2/tau^2;
    VarS = alpha*(1-alpha)/N; Mmax = 1/N;
    Bern(nn) = 2*exp(-(tau*alpha)^2/2/(VarS + Mmax*tau*alpha/3));
end
semilogy(Nlist, max(PGc,1e-300), [cMC 'o-'], 'LineWidth', lw, 'MarkerSize', ms), hold on
semilogy(Nlist, Cheb, [c00 '--'], 'LineWidth', lw)
semilogy(Nlist, Bern, [c01 '-.'], 'LineWidth', lw)
ylim([1e-20 1]), xlabel('N'), ylabel('P(G^c)')
legend('exact P(G^c)','Chebyshev \rho^2/\tau^2','Bernstein')
panel_label('(b)'), apply_fig_style(gca)
