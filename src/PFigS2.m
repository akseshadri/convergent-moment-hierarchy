%% PFigS2.m: Exact verification of the mixed-moment formulas (Supplementary Figure S2)
%
% A small heterogeneous, correlated, nonuniformly weighted network (N = 9).
% The closed forms of mu12, mu21, and kappa are compared with exact
% averages over all 2^9 reporting patterns.

clear, close all

N = 9; beta = (1:N)'/sum(1:N);
mi = 0.6 + 0.08*(1:N)' + 0.15*sin((1:N)');
sd = linspace(0.7,1.3,N)';
Sigma = diag(sd)*(0.45.^abs((1:N)'-(1:N)))*diag(sd);
alist = 0.1:0.1:1; % alpha
cf = zeros(numel(alist),3); ex = cf;
for ia = 1:numel(alist)
    cf(ia,:) = mixedmom(beta,mi,Sigma,alist(ia));
    [~, ~, ex(ia,:)] = enumvar(beta,mi,Sigma,alist(ia),[0 1]);
end
fprintf('max discrepancy relative to scale of each moment: %.2e\n', max(max(abs(cf-ex))./max(abs(ex))));

figure(102), set(gcf,'PaperUnits','inches','PaperPosition',[0 0 6.5 9])
lab = {'\mu_{1,2}','\mu_{2,1}','\kappa'};
for j = 1:3
    subplot(3,1,j), plot(alist,cf(:,j),'b-'), hold on, plot(alist,ex(:,j),'kx')
    ylabel(lab{j}), title(['(' char(96+j) ')'])
end
xlabel('\alpha'), subplot(3,1,1), legend('closed form','exact enumeration','Location','northeast'), legend boxoff
print('-depsc2','../figures/FigS02_mixed_moments_exact.eps')
