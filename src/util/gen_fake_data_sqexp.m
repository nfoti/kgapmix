% Generate fake data similar to the sngp paper for comparison.
% This code is based on similar code from Vinayak Rao's sngp code.

num_times = 10;
times = 0.1.*(1:num_times)';
kern_tau = 1/.1;
mus = times(2:2:num_times); % Atoms live at every other time point

nclusts = 12;

means = randn(1,nclusts)*10;
means(4) = -.2; % These are from Vinayak's code for some identifiability
means(5) = .3;
phis = ones(1,nclusts); % Could allow variable precisions for each cluster
Pi = gamrnd(1,1, 1,nclusts) + eps;
mu_inds = randsample(numel(mus), nclusts, 1);

% Number of observations per time point
nobs = 50;

Y = zeros(nobs*num_times,1);
X = zeros(nobs*num_times,1);
Z = zeros(nobs*num_times,1);

KK = sekern(times, mus, mu_inds, kern_tau, ones(nclusts,1));

for t = 1:num_times

  inds = (t-1)*nobs + (1:nobs);

  KKpi = KK(t,:).*Pi;
  prob = bsxfun(@rdivide, KKpi, sum(KKpi));

  Z(inds) = discreternd(nobs, prob);
  Y(inds) = rand(nobs,1)./phis(Z(inds))' + means(Z(inds))';

  X(inds) = times(t).*ones(nobs,1);

end
