function lp = logP(Y,KK,params)
% LOGP Compute log density of model.
%
% lp = logP()
%

N = numel(Y);

S = params.S;
Pi = params.Pi;
theta = params.theta;
phi = params.phi;
mus = params.mus;
psis = params.psis;
mu_inds = params.mu_inds;
psi_inds = params.psi_inds;

P = size(mus,1);
Q = numel(psis);

Ki = unique(S);
K = numel(Ki);

a0 = params.a0;
b0 = params.b0;
c0 = params.c0;
d0 = params.d0;
thtau0 = params.thtau0;
thmean0 = params.thmean0;
U = params.U;
V = params.V;

lp = 0;

% pi
lPi = log(Pi(Ki) + eps);
lp = lp + K*(a0*log(b0) - gammaln(a0)) + ...
          sum( (a0-1)*lPi - b0*Pi(Ki) );

%{
% These are wrong for mu and psi, should be like Umu(mu_inds(Ki)
% mu
if numel(Umu) == 1
  lp = lp + P*log(Umu+eps);
else
  lp = lp + sum(log(Umu+eps));
end

% psi
if numel(Upsi) == 1
  lp = lp + Q*log(Upsi+eps);
else
  lp = lp = sum(log(Upsi+eps));
end
%}

% theta
lp = lp - K*( 0.5*log(2*pi) + log(sqrt(thtau0)) ) ...
     - sum( 0.5*thtau0.*(theta(Ki) - thmean0).^2 );

% phi
lphi = log(phi(Ki) + eps);
lp = lp + K*(c0*log(d0) - gammaln(c0)) + ...
          sum( (c0-1)*lphi - d0*phi(Ki) );

% s
tmp = bsxfun(@times, KK, Pi);
tmp = bsxfun(@rdivide, tmp, sum(tmp,2));
inds = sub2ind(size(KK),(1:N)',S);
lp = lp + sum( log(tmp(inds)+eps) );

% data
lp = lp - N*0.5*log(2*pi) ...
        + sum( log(sqrt(phi(S))) - ...
               0.5*phi(S).*(Y' - theta(S)).^2 );

