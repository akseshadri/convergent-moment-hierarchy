%% PFig2.m: Validation with the IMD rainfall field (Figure 2)
%
% Uses the site means and covariance of the 1901-2011 IMD daily rainfall
% record as population moments, with uniform weights.
%   (a) exact variance and the (0,0), (0,1), (1,0), (1,1) approximations
%   (b) their relative errors with respect to the total variance
%   (c) their errors relative to the contribution of missing data, which is
%       Var(rhat) - Var(rbar), where rbar is the average with complete reporting

clear, close all

[mi, Sigma] = loadimd('../data/indiadat.mat');
N = numel(mi);
Vbar = sum(Sigma(:))/N^2;                 % variance of the complete-reporting average
alist = 0.1:0.1:1; % alpha
ord = [0 0; 0 1; 1 0; 1 1];
v = zeros(numel(alist),1); mu = zeros(numel(alist),4);
for ia = 1:numel(alist)
    v(ia) = varexact(mi,Sigma,alist(ia));
    for j = 1:4
        mu(ia,j) = truncvar(mi,Sigma,alist(ia),ord(j,1),ord(j,2));
    end
end
err = abs(mu - repmat(v,1,4));
ia = alist < 1;                           % all approximations are exact at alpha = 1
fprintf('alpha = 0.1: relative errors %.4f %.4f %.4f %.4f\n', err(1,:)/v(1));
fprintf('alpha = 0.1: Var = %.2f, Var(rbar) = %.2f, missing-data contribution = %.2f\n', v(1), Vbar, v(1)-Vbar);
Es2 = trace(Sigma)/N + mean((mi-mean(mi)).^2) - Vbar;   % expected spatial variance of a realization
fprintf('lambda = Var(rbar)/E s^2 = %.3f\n', Vbar/Es2);

%% Figure 2
figure(2), set(gcf,'PaperUnits','inches','PaperPosition',[0 0 13 4.2])
sty = {'r+-','bo-','md-','cs-'};
subplot(1,3,1), plot(alist,v,'k-','LineWidth',2), hold on
for j = 1:4, plot(alist,mu(:,j),sty{j}), end
xlabel('\alpha'), ylabel('variance (mm/day)^2'), title('(a)')
legend('exact','(0,0)','(0,1)','(1,0)','(1,1)','Location','northeast'), legend boxoff
subplot(1,3,2)
for j = 1:4, semilogy(alist(ia),err(ia,j)./v(ia),sty{j}), hold on, end
xlabel('\alpha'), ylabel('error relative to total variance'), grid on, title('(b)')
subplot(1,3,3)
for j = 1:4, semilogy(alist(ia),err(ia,j)./(v(ia)-Vbar),sty{j}), hold on, end
xlabel('\alpha'), ylabel('error relative to missing-data contribution'), grid on, title('(c)')
print('-depsc2','../figures/Fig02_imd_validation.eps')
