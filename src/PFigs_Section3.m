%% PFigs_Section3.m — Section 3 figures: closed forms and their limits
%
% Companion code for:
%   Seshadri, A. K. and Pal Majumder, A. (2026). A convergent moment
%   hierarchy for the variance of weighted spatial-mean estimators under
%   missing-at-random sampling.
%
% Validates the Section-3 closed-form variance truncations on the IMD rainfall
% field and probes the regimes where low-order truncations begin to fail.
%
%   Figure 3 (mc vs truncations): Monte Carlo vs the (0,0),(0,1),(1,0),(1,1)
%     truncations on the raw field, across the availability sweep alphalist.
%   Figure 4 (mixed moments): the mixed moments mu_{1,2}, mu_{2,1}, kappa,
%     Monte Carlo vs theory.
%   Figure 5 (decomposition, amplified means): term-by-term decomposition of
%     each truncation for a mean-amplified field (c = 20), where (0,0) visibly
%     separates from the others.
%   Figure 6 (stress regimes):
%     (a) long-range spatial dependence: Var(rhat) stops decaying with N;
%     (b) small networks: signed relative error of (0,1),(1,1),(0,2) against an
%         EXACT true variance (no Monte Carlo), on a symmetric-log axis, so the
%         genuine sign change of the (0,2) error is shown honestly rather than
%         as a spurious dip in its magnitude;
%     (c) concentrated weights: rho inflates and (1,1) recovers accuracy.
%
% Uses: truncations_scenarioI, truncations_general, ratio_var_exact_uniform,
% mc_truevar, apply_fig_style, panel_label.
%
% Data: 1 x 1 deg IMD gridded daily rainfall (../data/indiadat.mat).
%
% Requirements: Statistics and Machine Learning Toolbox (binornd); MATLAB
% R2016b or later (yyaxis); Parallel Computing Toolbox optional.


clear                                   

lw = 1.5; ms = 8;
cMC = 'k'; c00 = 'r'; c01 = 'b'; c10 = 'm'; c11 = 'c'; c02 = [0 0.6 0];

%% Load data
inddat = load('../data/indiadat.mat');
datmat = NaN(111,357,365);
for i = 1901:2011
    datmat(i-1900,:,:) = single(cell2mat(inddat.indiarainmodel(i)));
end
data0 = reshape(datmat,[357 111*365])';  

N = size(data0,2);
beta = ones(N,1)*1/N;                      % uniform weights (Scenario I)
alphalist = 0.1:0.1:1; Nalpha = numel(alphalist);   % reporting-probability (availability) sweep
Kensemb = 100000; Tuse = 365;              % Monte-Carlo ensemble size; days used (one year)

%% ---- helper anonymous transform (mean-anomaly, factor c and scale lambda) ----
meandatai = mean(data0,1); meandatamati = repmat(meandatai,[size(data0,1) 1]); mbar = mean(meandatai);
transform = @(lambda,c) mbar + c*(meandatamati-mbar) + lambda*(data0-meandatamati);

%% ================= FIGURES 3 & 4 : validation on the raw field =================
datause = data0(1:Tuse,:);
mi = mean(datause,1)'; sigmai2 = var(datause,[],1)'*(Tuse-1)/Tuse; sigmamat = cov(datause)*(Tuse-1)/Tuse;
m_r = mean(mi); Var_m = mean((mi-m_r).^2);
Sdiag = sum(sigmai2); Soff = sum(sum(sigmamat))-sum(diag(sigmamat)); M2m = sum(mi.^2);

varmc = NaN(Nalpha,1); mu21mc = varmc; mu12mc = varmc; kappamc = varmc;
mu00 = varmc; mu01 = varmc; mu10 = varmc; mu11 = varmc; mu12 = varmc; mu21 = varmc; kappal = varmc;
for k = 1:Nalpha
    alpha = alphalist(k);
    [varmc(k), mu21mc(k), mu12mc(k), kappamc(k)] = mc_truevar(datause, beta, alpha, Kensemb);
    T = truncations_scenarioI(alpha, N, Sdiag, Soff, M2m, Var_m, m_r);
    mu00(k)=T.mu00; mu01(k)=T.mu01; mu10(k)=T.mu10; mu11(k)=T.mu11;
    mu12(k)=T.mu12; mu21(k)=T.mu21; kappal(k)=T.kappa;
end

% ---- Fig 3: MC vs the four truncations ----
figure(3); clf
plot(alphalist, varmc, [cMC 'x-'], 'LineWidth', 2, 'MarkerSize', ms), hold on
plot(alphalist, mu00, [c00 '+-'], 'LineWidth', 1, 'MarkerSize', ms)
plot(alphalist, mu01, [c01 'o-'], 'LineWidth', 1, 'MarkerSize', ms)
plot(alphalist, mu10, [c10 'd-'], 'LineWidth', 1, 'MarkerSize', ms)
plot(alphalist, mu11, [c11 'diamond-'], 'LineWidth', 1, 'MarkerSize', ms)
xlabel('\alpha'), ylabel('\mu_2')
legend('MC','\mu_2^{(0,0)}','\mu_2^{(0,1)}','\mu_2^{(1,0)}','\mu_2^{(1,1)}')
apply_fig_style(gca)

% ---- Fig 4: mixed moments (MC vs theory) ----
figure(4); clf
subplot(3,1,1), plot(alphalist, mu12mc, [cMC 'x-'], 'LineWidth', 2), hold on, plot(alphalist, mu12, [c01 'o-'], 'LineWidth', lw)
    xlabel('\alpha'), ylabel('\mu_{1,2}'), legend('MC','theory'), panel_label('(a)'), apply_fig_style(gca)
subplot(3,1,2), plot(alphalist, mu21mc, [cMC 'x-'], 'LineWidth', 2), hold on, plot(alphalist, mu21, [c01 'o-'], 'LineWidth', lw)
    xlabel('\alpha'), ylabel('\mu_{2,1}'), legend('MC','theory'), panel_label('(b)'), apply_fig_style(gca)
subplot(3,1,3), plot(alphalist, kappamc, [cMC 'x-'], 'LineWidth', 2), hold on, plot(alphalist, kappal, [c01 'o-'], 'LineWidth', lw)
    xlabel('\alpha'), ylabel('\kappa'), legend('MC','theory'), panel_label('(c)'), apply_fig_style(gca)

%% ================= FIGURE 5 : decomposition, amplified means =================
datause = transform(0.5,20); datause = datause(1:Tuse,:);
mi = mean(datause,1)'; sigmai2 = var(datause,[],1)'*(Tuse-1)/Tuse; sigmamat = cov(datause)*(Tuse-1)/Tuse;
m_r = mean(mi); Var_m = mean((mi-m_r).^2);
Sdiag = sum(sigmai2); Soff = sum(sum(sigmamat))-sum(diag(sigmamat)); M2m = sum(mi.^2);

varmc = NaN(Nalpha,1); mu00 = varmc; mu01 = varmc; mu10 = varmc; mu11 = varmc;
T00 = NaN(Nalpha,3); T01 = NaN(Nalpha,3); T10 = NaN(Nalpha,5); T11 = NaN(Nalpha,6);
for k = 1:Nalpha
    alpha = alphalist(k);
    varmc(k) = mc_truevar(datause, beta, alpha, Kensemb);
    T = truncations_scenarioI(alpha, N, Sdiag, Soff, M2m, Var_m, m_r);
    mu00(k)=T.mu00; mu01(k)=T.mu01; mu10(k)=T.mu10; mu11(k)=T.mu11;
    T00(k,:)=T.t00; T01(k,:)=T.t01; T10(k,:)=T.t10; T11(k,:)=T.t11;
end
allterms = [T00(:); T01(:); T10(:); T11(:)];
pad = 0.1*(max(allterms)-min(allterms)); ylims = [min(allterms)-pad, max(allterms)+pad];

figure(5); clf
subplot(3,2,1), histogram(mi,20), hold on, histogram(sqrt(sigmai2),20)
    xlabel('value'), ylabel('count'), legend('m_i','\sigma_i'), panel_label('(a)'), apply_fig_style(gca)
subplot(3,2,2), plot(alphalist, varmc, [cMC 'x-'], 'LineWidth', 2, 'MarkerSize', ms), hold on
    plot(alphalist, mu00, [c00 '+-']), plot(alphalist, mu01, [c01 'o-']), plot(alphalist, mu10, [c10 'd-']), plot(alphalist, mu11, [c11 'diamond-'])
    xlabel('\alpha'), ylabel('\mu_2'), legend('MC','(0,0)','(0,1)','(1,0)','(1,1)'), panel_label('(b)'), apply_fig_style(gca)
subplot(3,2,3), plot(alphalist, T00, 'LineWidth', lw), xlabel('\alpha'), ylabel('\mu_2^{(0,0)} terms'), ylim(ylims)
    legend('1','2','3'), panel_label('(c)'), apply_fig_style(gca)
subplot(3,2,4), plot(alphalist, T01, 'LineWidth', lw), xlabel('\alpha'), ylabel('\mu_2^{(0,1)} terms'), ylim(ylims)
    legend('1','2','3'), panel_label('(d)'), apply_fig_style(gca)
subplot(3,2,5), plot(alphalist, T10, 'LineWidth', lw), xlabel('\alpha'), ylabel('\mu_2^{(1,0)} terms'), ylim(ylims)
    legend('1','2','3','4','5'), panel_label('(e)'), apply_fig_style(gca)
subplot(3,2,6), plot(alphalist, T11, 'LineWidth', lw), xlabel('\alpha'), ylabel('\mu_2^{(1,1)} terms'), ylim(ylims)
    legend('1','2','3','4','5','6'), panel_label('(f)'), apply_fig_style(gca)

%% ================= FIGURE 6 : stress regimes =================
sigma0 = sqrt(mean(var(data0(1:Tuse,:),[],1)*(Tuse-1)/Tuse));
Kstress = 8000; astress = 0.4;               % MC size and fixed availability for the stress test
figure(6); clf

% ---- (a) long-range dependence: variance vs N stops decaying ----
subplot(1,3,1)
Nlist = [20 40 80 160 300];
rclist = [0 0.85];                           % shared-mode fraction: 0 = independent sites, 0.85 = strong long-range dependence
rcsty = {[cMC 'o-'], [cMC 's-']};
for rr = 1:numel(rclist)
    rc = rclist(rr); vv = NaN(numel(Nlist),1);
    for nn = 1:numel(Nlist)
        N = Nlist(nn); bU = ones(N,1)/N; mi = (sigma0*linspace(0.6,1.4,N))';
        g = randn(Tuse,1); e = randn(Tuse,N);
        fld = repmat(mi',Tuse,1) + sigma0*(sqrt(rc)*repmat(g,1,N) + sqrt(1-rc)*e);
        vv(nn) = mc_truevar(fld, bU, astress, Kstress);
    end
    loglog(Nlist, vv, rcsty{rr}, 'LineWidth', lw, 'MarkerSize', 6), hold on
end
xlabel('N'), ylabel('Var(rhat)'), legend('weak dependence','strong long-range')
panel_label('(a)'), apply_fig_style(gca)

% ---- (b) small network: (0,1) degrades, (0,2) improves ----
% Fixed population field (deterministic moments; independent sites) with an
% EXACT true variance (closed form, no Monte-Carlo noise), so each curve is a
% smooth function of N.  The (0,2) truncation error changes sign with N, so we
% plot the signed relative error on a symmetric-log axis
subplot(1,3,2)
Nlist = [200 120 80 50 30 18 10 6];
s01 = NaN(numel(Nlist),1); s11 = s01; s02 = s01;
for nn = 1:numel(Nlist)
    N  = Nlist(nn);
    bU = ones(N,1)/N;
    mi = (sigma0*linspace(0.6,1.4,N))';    % fixed heterogeneous mean profile
    s2 = sigma0^2*ones(N,1);               % common site variance (independent sites)
    T  = truncations_general(bU, mi, diag(s2), astress);   % exact population moments
    v  = ratio_var_exact_uniform(mi, s2, astress);         % exact variance, no Monte Carlo
    s01(nn) = (v - T.mu01)/v;
    s11(nn) = (v - T.mu11)/v;
    s02(nn) = (v - T.mu02)/v;
end
Csl = 1e-3;                                  % symmetric-log linear threshold
sy  = @(y) sign(y).*log10(1 + abs(y)/Csl);   % symlog: linear near 0, log in the tails
plot([max(Nlist) min(Nlist)], [0 0], ':', 'Color', [.6 .6 .6], 'HandleVisibility','off'), hold on
plot(Nlist, sy(s01), [c01 'o-'], 'LineWidth', lw)
plot(Nlist, sy(s11), [c11 's-'], 'LineWidth', lw)
plot(Nlist, sy(s02), 'Color', c02, 'Marker', '^', 'LineStyle', '-', 'LineWidth', lw)
yt = [-1e-1 -1e-2 -1e-3 0 1e-3 1e-2 1e-1];
set(gca, 'YTick', sy(yt), 'YTickLabel', ...
    {'-10^{-1}','-10^{-2}','-10^{-3}','0','10^{-3}','10^{-2}','10^{-1}'})
set(gca, 'XScale','log', 'XDir','reverse')
xlabel('network size N (harder \rightarrow)'), ylabel('signed relative error')
legend('(0,1)','(1,1)','(0,2)','Location','southwest')
panel_label('(b)'), apply_fig_style(gca)

% ---- (c) concentrated weights: rho inflates, (1,1) improves ----
subplot(1,3,3)
N = 60; mi = (sigma0*linspace(0.6,1.4,N))';
fld = repmat(mi',Tuse,1) + sigma0*randn(Tuse,N); sigmamat = cov(fld)*(Tuse-1)/Tuse; mihat = mean(fld,1)';
wlist = [1/N 0.05 0.1 0.2 0.35 0.5];         % weight on site 1: 1/N (uniform) up to 0.5 (concentrated)
e01 = NaN(numel(wlist),1); e11 = e01; e02 = e01; rho_list = e01;
for ww = 1:numel(wlist)
    w = wlist(ww); bW = NaN(N,1); bW(1)=w; bW(2:N)=(1-w)/(N-1);
    T = truncations_general(bW, mihat, sigmamat, astress);
    v = mc_truevar(fld, bW, astress, Kstress);
    e01(ww)=abs(v-T.mu01)/abs(v); e11(ww)=abs(v-T.mu11)/abs(v); e02(ww)=abs(v-T.mu02)/abs(v); rho_list(ww)=T.rho;
end
semilogy(wlist, e01, [c01 'o-'], 'LineWidth', lw), hold on
semilogy(wlist, e11, [c11 's-'], 'LineWidth', lw)
semilogy(wlist, e02, 'Color', c02, 'Marker', '^', 'LineStyle', '-', 'LineWidth', lw)
yyaxis right, plot(wlist, rho_list, [c00 ':'], 'LineWidth', lw), ylabel('\rho')
yyaxis left, xlabel('weight on site 1, \beta_1 (uniform \rightarrow concentrated)'), ylabel('relative error')
legend('(0,1)','(1,1)','(0,2)','\rho'), panel_label('(c)'), apply_fig_style(gca)
