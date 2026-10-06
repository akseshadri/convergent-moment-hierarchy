%% PFig1.m: Error scaling and denominator concentration (Figure 1)
%
% This figure uses the homogeneous benchmark: i.i.d. sites with m = sigma = 1 and uniform weights.
%   (a) exact errors of the (0,0) and (0,1) approximations against rho, with
%       reference slopes 2 and 4 (fixed alpha = 0.3, N varies)
%   (b) relative errors against N for alpha = 0.1, 0.5, 0.9
%   (c) distribution of the normalized denominator (S - ES)/ES, with the good
%       event |S - ES| <= tau ES shaded (alpha = 0.5, tau = 0.5)
%   (d) probability of the rare event G^c, with Chebyshev and Bernstein bounds

clear, close all

%% (a) error against rho at fixed alpha
alpha = 0.3;
Nlist = [50 80 120 200 320 500 800 1200 1800];
rho = sqrt((1-alpha)./(alpha*Nlist));
e00 = zeros(size(Nlist)); e01 = e00;
for i = 1:numel(Nlist)
    N = Nlist(i); mi = ones(N,1); Sigma = eye(N);
    v = varexact(mi,Sigma,alpha);
    e00(i) = abs(v - truncvar(mi,Sigma,alpha,0,0));
    e01(i) = abs(v - truncvar(mi,Sigma,alpha,0,1));
end

%% (b) relative error against N
alist = [0.1 0.5 0.9]; % alpha
Nlist2 = [50 80 120 200 400 800];
r00 = zeros(numel(alist),numel(Nlist2)); r01 = r00;
for ia = 1:numel(alist)
    for i = 1:numel(Nlist2)
        N = Nlist2(i); mi = ones(N,1); Sigma = eye(N);
        v = varexact(mi,Sigma,alist(ia));
        r00(ia,i) = abs(v - truncvar(mi,Sigma,alist(ia),0,0))/v;
        r01(ia,i) = abs(v - truncvar(mi,Sigma,alist(ia),0,1))/v;
    end
end

%% (c)-(d) denominator concentration (S = K/N under uniform weights)
alpha = 0.5; tau = 0.5;
Nlist3 = [25 50 100 200 357 600 1000];
PG = zeros(size(Nlist3)); cheb = PG; bern = PG;
for i = 1:numel(Nlist3)
    N = Nlist3(i); k = (0:N)';
    pk = exp(gammaln(N+1)-gammaln(k+1)-gammaln(N-k+1)+k*log(alpha)+(N-k)*log(1-alpha));
    PG(i) = sum(pk(abs(k/N-alpha) > tau*alpha));
    cheb(i) = min(1, (1-alpha)/(alpha*N)/tau^2);
    bern(i) = min(1, 2*exp(-(tau*alpha)^2/2/(alpha*(1-alpha)/N + tau*alpha/(3*N))));
end

%% Figure 1
figure(1), set(gcf,'PaperUnits','inches','PaperPosition',[0 0 10 8.5])
subplot(2,2,1), loglog(rho,e00,'ro-','LineWidth',1.0), hold on, loglog(rho,e01,'bs-','LineWidth',1.0)
loglog(rho,rho.^2,'r--'), loglog(rho,e01(end)*(rho/rho(end)).^4,'b--')
xlabel('\rho (via N, \alpha = 0.3)'), ylabel('exact absolute variance error')
legend('(0,0)','(0,1)','slope 2','slope 4','Location','southeast'), legend boxoff, grid on, title('(a)')
subplot(2,2,2), sty = {'o','s','d'};
for ia = 1:numel(alist)
    loglog(Nlist2,r00(ia,:),['k-' sty{ia}]), hold on
    loglog(Nlist2,r01(ia,:),['b--' sty{ia}])
end
xlabel('N'), ylabel('exact relative error'), grid on, title('(b)')
legend('\alpha=.1 (0,0)','\alpha=.1 (0,1)','\alpha=.5 (0,0)','\alpha=.5 (0,1)','\alpha=.9 (0,0)','\alpha=.9 (0,1)','Location','southwest'), legend boxoff
subplot(2,2,3), lsty = {'b-','r-','k-'}; Nc = [50 200 357];
patch([-tau tau tau -tau],[0 0 8 8],[0.9 0.9 0.9],'EdgeColor','none'), hold on
for i = 1:numel(Nc)
    N = Nc(i); k = (0:N)';
    pk = exp(gammaln(N+1)-gammaln(k+1)-gammaln(N-k+1)+k*log(alpha)+(N-k)*log(1-alpha));
    plot((k/N-alpha)/alpha, pk*N*alpha, lsty{i}, 'LineWidth',1.0)
end
xlim([-1 1]), ylim([0 8]), xlabel('(S - E S)/E S'), ylabel('exact discrete density'), title('(c)')
legend('good event','N=50','N=200','N=357','Location','northeast'), legend boxoff
subplot(2,2,4), semilogy(Nlist3,max(PG,realmin),'ko-','LineWidth',1.0), hold on
semilogy(Nlist3,cheb,'r--'), semilogy(Nlist3,bern,'b-.')
ylim([1e-60 1]), set(gca,'YTick',10.^(-60:10:0)), xlabel('N'), ylabel('P(G^c)'), grid on, title('(d)')
legend('exact','Chebyshev','Bernstein','Location','southwest'), legend boxoff
print('-depsc2','../figures/Fig01_scaling_concentration.eps')
