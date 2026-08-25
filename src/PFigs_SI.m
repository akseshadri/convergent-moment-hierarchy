%% PFigs_SI.m — Supplementary figures: variance decomposition
%
% Companion code for:
%   Seshadri, A. K. and Pal Majumder, A. (2026). A convergent moment
%   hierarchy for the variance of weighted spatial-mean estimators under
%   missing-at-random sampling.
%
% Term-by-term decomposition of the variance truncations in the two regimes
% not shown in the main text: the raw rainfall field (site means
% negligible relative to variability, so all truncations coincide) and a
% reduced-temporal-variance case (the (0,0) and (0,1) truncations first begin
% to separate). The decisive mean-amplified case is main-text Figure 5.
%
%   Figure S1 (decomposition, raw):     lambda = 1,   c = 1.
%   Figure S2 (decomposition, reduced): lambda = 0.5, c = 1.
%
% The field transform is
%   x -> mbar + c*(site_mean - mbar) + lambda*(x - site_mean),
% where c scales the between-site mean contrast and lambda scales the temporal
% (within-site) variability.
%
% Uses: truncations_scenarioI, mc_truevar, apply_fig_style, panel_label.
%
% Data: 1 x 1 deg IMD gridded daily rainfall (../data/indiadat.mat).
%
% Requirements: Statistics and Machine Learning Toolbox (binornd); Parallel
% Computing Toolbox optional.


clear                                   % keep any open figures; run_all_figures closes them once at the start

lw = 1.5; ms = 8;
cMC = 'k'; c00 = 'r'; c01 = 'b'; c10 = 'm'; c11 = 'c';

%% Load data
inddat = load('../data/indiadat.mat');
datmat = NaN(111,357,365);
for i = 1901:2011
    datmat(i-1900,:,:) = single(cell2mat(inddat.indiarainmodel(i)));
end
data0 = reshape(datmat,[357 111*365])';   % [time,location], raw field
% (portable alternative:  data0 = readmatrix('data.txt');)

N = size(data0,2);
beta = ones(N,1)*1/N;
alphalist = 0.1:0.1:1; Nalpha = numel(alphalist);
Kensemb = 100000; Tuse = 365;             % Monte-Carlo ensemble size; days used (one year)

meandatai = mean(data0,1); meandatamati = repmat(meandatai,[size(data0,1) 1]); mbar = mean(meandatai);
transform = @(lambda,c) mbar + c*(meandatamati-mbar) + lambda*(data0-meandatamati);

%% ---- Supplementary settings: (lambda, c, SI number, short name) ----
lamlist = [1   0.5];   % temporal-variability scale lambda for S1, S2
clist   = [1   1  ];   % between-site mean-contrast scale c for S1, S2

for cc = 1:numel(clist)

    lambda = lamlist(cc); c = clist(cc);
    datause = transform(lambda,c); datause = datause(1:Tuse,:);
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

    figure(6+cc); clf   % S1 -> figure(7), S2 -> figure(8): distinct numbers, cleared each pass
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

end
