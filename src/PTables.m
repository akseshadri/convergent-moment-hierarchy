%% PTables.m: Values for Table 4 and Supplementary Table S3
%
% Homogeneous benchmark: i.i.d. sites with m = sigma = 1 and uniform weights.
% Table 4: convergence at fixed N = 10.
% Table S3: errors along the shift-invariant truncations (n, n+1) at alpha = 0.3.

clear, close all

N = 10; mi = ones(N,1); Sigma = eye(N);
fprintf('Table 4(a): mu^(n,n+1), N = %d\n', N);
nlist = [5 10 20 40 41 80 81];
for alpha = [0.45 0.50 0.55]
    fprintf('alpha = %.2f  exact = %.4f |', alpha, varexact(mi,Sigma,alpha));
    for n = nlist, fprintf(' %.4g', truncvar(mi,Sigma,alpha,n,n+1)); end
    fprintf('\n');
end
alpha = 0.6; nlist = [1 2 4 8 16 32 64];
fprintf('Table 4(b): alpha = 0.6, exact = %.4f\n  (0,n):  ', varexact(mi,Sigma,alpha));
for n = nlist, fprintf(' %.4f', truncvar(mi,Sigma,alpha,0,n)); end
fprintf('\n  (n,n+1):');
for n = nlist, fprintf(' %.4f', truncvar(mi,Sigma,alpha,n,n+1)); end
fprintf('\n');

fprintf('Table S3: |Var - mu^(n,n+1)| at alpha = 0.3\n');
alpha = 0.3;
for N = [40 80 160 320 640]
    mi = ones(N,1); Sigma = eye(N); v = varexact(mi,Sigma,alpha);
    fprintf('N = %4d', N);
    for n = 0:4, fprintf('  %.3e', abs(truncvar(mi,Sigma,alpha,n,n+1) - v)); end
    fprintf('\n');
end
