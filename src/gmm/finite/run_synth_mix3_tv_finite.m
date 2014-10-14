% Run truncated gibbs sampler on toy mixture of 3 components that vary with
% time

s = RandStream('mt19937ar', 'Seed', 8675309);
%s = RandStream('mt19937ar', 'Seed', floor(rand*1000000));
RandStream.setGlobalStream(s);
fprintf('RandStream Seed: %d\n', RandStream.getGlobalStream.Seed);

% Generate data
N = 250;
Y = [];
XX = [];
Z = [];
means = [];
taus = [];
p = [];
gen_mix3_tv;

spread = max(Y) - min(Y);
mY = mean(Y);
vbar = var(Y);

KK = 10;
a0 = 1/KK;
b0 = 1;

% thmean0 = data mean
% thtau0  = .5/sample variance
% c0 and d0 are chosen so that mean of phi is 10/vbar and variance is
%  50/vbar (taken from Yee Whye's MCMC for NRM paper)
thmean0 = mY;
thtau0 = .5/vbar;
c0 = 2;
d0 = vbar/5;
U = 1;
V = 1;

S = randsample(KK,N,true);
Pi = a0/b0*ones(1,KK);
theta = randn(1,KK)./sqrt(thtau0) + thmean0;
phi = gamrnd(c0, 1/d0, 1, KK);

mus = unique(XX);
psis = 1;
mu_inds = randsample(size(mus,1),KK,true);
psi_inds = ones(1,KK);

nburn = 1;
thin = 3000;

init_burn = init_params_struct(Y, XX, 'nsamp', nburn, 'thin', thin, ...
                               'S', S, 'Pi', Pi, 'theta', theta, ...
                               'phi', phi, 'mus', mus, 'psis', psis, ...
                               'mu_inds', mu_inds, 'psi_inds', psi_inds, ...
                               'thmean0', thmean0, 'thtau0', thtau0, ...
                               'a0', a0, 'b0', b0, 'c0', c0, 'd0', d0, ...
                               'U', U, 'V', V, 'kernfun', @sekern ...
                              );
fprintf('Burn in...\n');                   
params = kgap_gmm_finite(Y,XX,init_burn);

lp_burn = params.lp;

nsamp = 100;
thin = 5;
init = init_params_struct(Y, XX, 'nsamp', nsamp, 'thin', thin, ...
                          'S', params.S, 'Pi', params.Pi, ...
                          'theta', params.theta, 'phi', params.phi, ...
                          'mus', mus, 'psis', psis, ...
                          'mu_inds', params.mu_inds, ...
                          'psi_inds', params.psi_inds, ...
                          'thmean0', thmean0, 'thtau0', thtau0, ...
                          'a0', a0, 'b0', b0, 'c0', d0, 'd0', d0, ...
                          'U', U, 'V', V, 'kernfun', @sekern ...
                         );

fprintf('Sampling...\n');
params = kgap_gmm_finite(Y,XX,init);
                       
% Evaluate
M = min(nsamp,50);
nc = zeros(1,M);
ari = zeros(1,M);
for m = 1:M
  nc(m) = numel(unique(params.S_samp(:,m)));
  ari(m) = adjrand(params.S_samp(:,m),Z);
end

fprintf('Num clusters: %.2f (+- %.2f)\n', mean(nc), std(nc));
fprintf('ARI with truth: %.3f (+- %.3f)\n', mean(ari), std(ari));


% Draw samples for distribution at covariate values as well as at
% interpolating values for the last sample (could average this over
% posterior draws as well)
Xev = linspace(1,3,4)';
Nev = 200;
ind = params.nsamp; % Just use last sample for now.
Ymin = min(Y);
Ymax = max(Y);
figure;
nrow = ceil(sqrt(size(Xev,1)));
for i = 1:size(Xev,1)
  subplot(nrow,nrow,i);
  KKpi = sekern(Xev(i),mus,params.mu_inds_samp(:,ind),psis,params.psi_inds_samp(:,ind)).*params.Pi_samp(:,ind)';
  KKpi = bsxfun(@rdivide, KKpi, sum(KKpi));
  xxx = linspace(Ymin,Ymax,500)';
  Fx = exp(-repmat(params.phi,N,1).*(repmat(xxx,1,KK)-repmat(params.theta,N,1)).^2) * KKpi';
  plot(xxx,Fx);
  ss = sprintf('X: %2.2f', Xev(i));
  title(ss);
end