%% PFigS3.m: Term-by-term decomposition for the IMD field (Supplementary Figure S3)
%
% Uniform weights.  (a) site means and standard deviations; (b) exact variance
% and the four low-order approximations; (c)-(f) the terms of (0,0), (0,1),
% (1,0), and (1,1), in the order of Supplementary Table S2.

clear, close all

[mi, Sigma] = loadimd('../data/indiadat.mat');
N = numel(mi); beta = ones(N,1)/N;
Sd = trace(Sigma); So = sum(Sigma(:)) - Sd;
Sm = sum(mi.^2); Vm = mean((mi-mean(mi)).^2); mr = mean(mi);
alist = 0.1:0.1:1; % alpha
nA = numel(alist);
v = zeros(nA,1); T00 = zeros(nA,3); T01 = T00; T10 = zeros(nA,5); T11 = zeros(nA,6);
for ia = 1:nA
    a = alist(ia); v(ia) = varexact(mi,Sigma,a);
    mm = mixedmom(beta,mi,Sigma,a);       % [mu12 mu21 kappa]
    t1 = Sd/(a*N^2); t2 = So/N^2; t3 = (1-a)/a*Sm/N^2; t3b = (1-a)/a*Vm/N;
    c4 = -2*mm(2)/a^3; c5 = mm(3)/a^4; c5b = 2*a*mr*mm(1)/a^4;
    T00(ia,:) = [t1 t2 t3]; T01(ia,:) = [t1 t2 t3b];
    T10(ia,:) = [t1 t2 t3 c4 c5]; T11(ia,:) = [t1 t2 t3b c4 c5b c5];
end

figure(103), set(gcf,'PaperUnits','inches','PaperPosition',[0 0 10 11])
subplot(3,2,1), plot(mi,sqrt(diag(Sigma)),'k.'), xlabel('site mean m_i'), ylabel('site SD \sigma_i'), title('(a)')
subplot(3,2,2), plot(alist,v,'k-','LineWidth',2), hold on
plot(alist,sum(T00,2),'r+-'), plot(alist,sum(T01,2),'bo-'), plot(alist,sum(T10,2),'md-'), plot(alist,sum(T11,2),'cs-')
xlabel('\alpha'), ylabel('variance'), title('(b)'), legend('exact','(0,0)','(0,1)','(1,0)','(1,1)','Location','northeast'), legend boxoff
subplot(3,2,3), plot(alist,T00), xlabel('\alpha'), ylabel('(0,0) terms'), title('(c)')
legend('site variances','spatial covariance','mean term','Location','northeast'), legend boxoff
subplot(3,2,4), plot(alist,T01), xlabel('\alpha'), ylabel('(0,1) terms'), title('(d)')
legend('site variances','spatial covariance','spread of site means','Location','northeast'), legend boxoff
subplot(3,2,5), plot(alist,T10), xlabel('\alpha'), ylabel('(1,0) terms'), title('(e)')
legend('site variances','spatial covariance','mean term','-2\mu_{2,1}/\alpha^3','\kappa/\alpha^4','Location','northeast'), legend boxoff
subplot(3,2,6), plot(alist,T11), xlabel('\alpha'), ylabel('(1,1) terms'), title('(f)')
legend('site variances','spatial covariance','spread of site means','-2\mu_{2,1}/\alpha^3','2m_R\mu_{1,2}/\alpha^4','\kappa/\alpha^4','Location','northeast'), legend boxoff
print('-depsc2','../figures/FigS03_decomposition_raw.eps')
