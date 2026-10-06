%% PFigS1.m: Homogeneous two-parameter scaling (Supplementary Figure S1)
%
% i.i.d. sites with m = sigma = 1 and uniform weights.
%   (a) exact (0,0) error against rho as both alpha and N vary, with the
%       leading reference m^2 rho^2
%   (b) |Var - mu^(0,1)|/rho^4 against N for several alpha; dashed lines
%       show the limit sigma^2/(1-alpha)

clear, close all

alist = [0.10:0.10:0.90]; %alpha

Nlist = [40 80 120 160 200 240 280 320];

figure(101), set(gcf,'PaperUnits','inches','PaperPosition',[0 0 10 4.2])
subplot(1,2,1)
for ia = 1:numel(alist)
    rho = sqrt((1-alist(ia))./(alist(ia)*Nlist)); e = zeros(size(Nlist));
    for i = 1:numel(Nlist)
        N = Nlist(i); v = varexact(ones(N,1),eye(N),alist(ia));
        e(i) = abs(v - truncvar(ones(N,1),eye(N),alist(ia),0,0));
    end
    loglog(rho,e,'o'), hold on
end
rr = logspace(-2,log10(0.5),20); loglog(rr,rr.^2,'k--')
xlabel('\rho'), ylabel('exact (0,0) error'), grid on, title('(a)')

alist = [0.15 0.30 0.50 0.80]; Nlist = [40 80 160 320 640];
subplot(1,2,2), sty = {'bo-','rs-','gd-','m^-'};
for ia = 1:numel(alist)
    a = alist(ia); q = zeros(size(Nlist));
    for i = 1:numel(Nlist)
        N = Nlist(i); v = varexact(ones(N,1),eye(N),a);
        q(i) = abs(v - truncvar(ones(N,1),eye(N),a,0,1))/((1-a)/(a*N))^2;
    end
    semilogx(Nlist,q,sty{ia}), hold on
    semilogx(Nlist,1/(1-a)*ones(size(Nlist)),[sty{ia}(1) '--'])
end
xlabel('N'), ylabel('|Var - \mu^{(0,1)}| / \rho^4'), grid on, title('(b)')
print('-depsc2','../figures/FigS01_homogeneous_scaling.eps')
