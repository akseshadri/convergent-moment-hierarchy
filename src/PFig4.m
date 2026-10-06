%% PFig4.m: Stress tests (Figure 4)
%
% Three ways in which the low-order approximations become less accurate,
% at fixed availability alpha = 0.4.
%   (a) Strong spatial correlation: uniform weights, N = 100, unit variances,
%       equal correlation c between all pairs of sites.  Error of (0,1)
%       relative to the total variance and to the contribution of missing
%       data, Var(rhat) - Var(rbar), with the leading-term prediction
%       rho^2 Var(rbar) of the exact decomposition.
%   (b) Small networks: uniform weights, independent unit-variance sites with
%       means spanning 0.6-1.4, N = 200 down to 6.  Relative errors against rho.
%   (c) Concentrated weights: N independent sites with unit means and
%       variances; site 1 carries weight w and the others share 1 - w.
%       Exact values computed by enumerating all 2^N reporting patterns.
%       Therefore N is chosen to be small e.g. N=20. 

clear, close all

alpha = 0.4;

%% (a) strong spatial correlation
N = 100; clist = 0.1:0.1:0.9;
rho2 = (1-alpha)/(alpha*N);
eTot = zeros(size(clist)); eMis = eTot; eLead = eTot; lam = eTot;
for i = 1:numel(clist)
    c = clist(i); Sigma = (1-c)*eye(N) + c*ones(N); mi = ones(N,1);
    v = varexact(mi,Sigma,alpha);
    Vbar = sum(Sigma(:))/N^2;
    e = truncvar(mi,Sigma,alpha,0,1) - v;
    eTot(i) = e/v; eMis(i) = e/(v - Vbar); eLead(i) = rho2*Vbar/(v - Vbar);
    lam(i) = Vbar/(trace(Sigma)/N - Vbar);  % diagnostic lambda = Var(rbar)/E s^2
end
fprintf('(a) c = 0.1, 0.5, 0.9: error/missing-data = %.3f %.3f %.3f, lambda = %.3f %.3f %.3f\n', eMis([1 5 9]), lam([1 5 9]));

%% (b) small networks, uniform weights
Nlist = [200:-10:20 18:-2:6 5];
ord = [0 1; 1 1; 0 2];
rhoB = sqrt((1-alpha)./(alpha*Nlist)); eB = zeros(numel(Nlist),3);
for i = 1:numel(Nlist)
    N = Nlist(i); mi = linspace(0.6,1.4,N)'; Sigma = eye(N);
    v = varexact(mi,Sigma,alpha);
    for j = 1:3
        eB(i,j) = abs(truncvar(mi,Sigma,alpha,ord(j,1),ord(j,2)) - v)/v;
    end
end

%% (c) concentrated weights, N = 20
N = 20; wlist = [0.004:0.004:0.02 0.04:0.02:0.10 0.20 0.30 0.40 0.50 0.60 0.70];
rhoC = zeros(size(wlist)); eC = zeros(numel(wlist),3);
for i = 1:numel(wlist)
    w = wlist(i); beta = [w; (1-w)/(N-1)*ones(N-1,1)];
    rhoC(i) = sqrt((1-alpha)/alpha*sum(beta.^2));
    [v, mu] = enumvar(beta, ones(N,1), eye(N), alpha, ord);
    eC(i,:) = abs(mu' - v)/v;
end
fprintf('(b) (0,1) errors: '); fprintf('%.3f ', eB(:,1)); fprintf('\n');
fprintf('(c) (0,1) errors: '); fprintf('%.3f ', eC(:,1)); fprintf('\n');

%% Figure 4
figure(4), set(gcf,'PaperUnits','inches','PaperPosition',[0 0 14 4.5])
subplot(1,3,1), semilogy(clist,abs(eTot),'bo-','LineWidth',1.0), hold on
semilogy(clist,abs(eMis),'rs-','LineWidth',1.0), semilogy(clist,eLead,'r--')
xlabel('correlation c'), ylabel('relative error of (0,1)'), grid on, title('(a)')
legend('relative to total','relative to missing-data part','leading-term (exact)','Location','northwest'), legend boxoff
sty = {'bo-','cs-','g^-'};
subplot(1,3,2)
for j = 1:3, semilogy(rhoB,eB(:,j),sty{j}), hold on, end
xlim([0.08 0.9]), ylim([1e-5 5]), xlabel('\rho'), ylabel('relative error'), grid on, title('(b) uniform weights, N = 200 to 5')
legend('(0,1)','(1,1)','(0,2)','Location','southeast'), legend boxoff
subplot(1,3,3)
for j = 1:3, semilogy(rhoC,eC(:,j),sty{j}), hold on, end
xlim([0.08 0.9]), ylim([1e-5 5]), xlabel('\rho'), ylabel('relative error'), grid on, title('(c) concentrated weights, N = 20')
print('-depsc2','../figures/Fig04_stress.eps')
