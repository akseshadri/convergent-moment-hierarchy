%% PFig3.m: Common shift of the field (Figure 3)
%
% A constant b is added to every site mean of the IMD field (alpha = 0.4,
% uniform weights).  The exact variance is unchanged by the shift.  Only the
% shift-invariant truncation (0,1) retains this property; variance of (0,0), (1,1), 
% and (0,2) spuriously depend on the mean level.
%   (a) exact variance and the four approximations against b/sigma_RMS
%   (b) their relative errors

clear, close all

[mi, Sigma] = loadimd('../data/indiadat.mat');
alpha = 0.4;
sigrms = sqrt(mean(diag(Sigma)));
shift = [0 0.5 1 2 4];
ord = [0 0; 0 1; 1 1; 0 2];
v = zeros(numel(shift),1); mu = zeros(numel(shift),4);
for i = 1:numel(shift)
    mib = mi + shift(i)*sigrms;
    v(i) = varexact(mib,Sigma,alpha);
    for j = 1:4
        mu(i,j) = truncvar(mib,Sigma,alpha,ord(j,1),ord(j,2));
    end
end

%% Figure 3
figure(3), set(gcf,'PaperUnits','inches','PaperPosition',[0 0 10 4.2])
sty = {'r+-','bo-','cs-','g^-'};
subplot(1,2,1), plot(shift,v,'k-','LineWidth',2), hold on
for j = 1:4, plot(shift,mu(:,j),sty{j}), end
xlabel('b/\sigma_{RMS}'), ylabel('variance (mm/day)^2'), title('(a)')
legend('exact','(0,0)','(0,1)','(1,1)','(0,2)','Location','northwest'), legend boxoff
subplot(1,2,2)
for j = 1:4, semilogy(shift,abs(mu(:,j)-v)./v,sty{j}), hold on, end
xlabel('b/\sigma_{RMS}'), ylabel('relative error'), grid on, title('(b)')
print('-depsc2','../figures/Fig03_translation.eps')
