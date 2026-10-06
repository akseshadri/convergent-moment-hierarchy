function [mi, Sigma, T] = loadimd(datafile)
% LOADIMD  Site means and covariance of the IMD daily rainfall field.
%
%   [mi, Sigma, T] = loadimd(datafile) loads the 1 x 1 deg gridded daily
%   rainfall over India for 1901-2011 (357 locations, 365 days per year),
%   stacks all days into a T x N record, and returns the site means mi
%   (N x 1), the covariance matrix Sigma (N x N, normalization 1/T), and the
%   number of days T.  These are treated as population moments of the field.

dat = load(datafile);
X = [];
for yr = 1901:2011
    Xi = dat.indiarainmodel{yr};          % 357 x 365 array (site x day)
    X = [X; Xi'];                         % append as day x site
end
T = size(X,1);
mi = mean(X,1)';
Xc = X - repmat(mi',T,1);
Sigma = (Xc'*Xc)/T;
