% Generate fake data similar to the sngp paper for comparison.
% This code is based on similar code from Vinayak Rao's sngp code.

num_times = 10;
times = 0.1.*(1:num_times)';
% Atoms active (on average for non-uniform times) for 2 time steps centered on current time
kern_tau = 2*mean(diff(times)); 
mus = times(2:2:num_times); % Atoms live at every other time point

nclusts = 12;

var_means = 10;
means = randn(1,nclusts)*var_means;
means(4) = -.2; % These are from Vinayak's code for some identifiability
means(5) = .3;
phis = ones(1,nclusts); % Could allow variable precisions for each cluster
Pi = gamrnd(1,1, 1,nclusts) + eps;
mu_inds = randsample(numel(mus), nclusts, 1);

% Number of observations per time point
nobs = 50;
nobs_te = 5;

Y = zeros(nobs*num_times,1);
X = zeros(nobs*num_times,1);
Z = zeros(nobs*num_times,1);

Ytest = zeros(nobs_te*num_times,1);
Xtest = zeros(nobs_te*num_times,1);
Ztest = zeros(nobs_te*num_times,1);

KK = boxkern(times, mus, mu_inds, kern_tau, ones(nclusts,1));

for t = 1:num_times

  inds = (t-1)*nobs + (1:nobs);
  inds_te = (t-1)*nobs_te + (1:nobs_te);

  KKpi = KK(t,:);%.*Pi; % Just uniformly sample from active atoms at each time
  prob = bsxfun(@rdivide, KKpi, sum(KKpi));

  Z(inds) = discreternd(nobs, prob);
  Y(inds) = randn(nobs,1)./phis(Z(inds))' + means(Z(inds))';

  X(inds) = times(t).*ones(nobs,1);

  Ztest(inds_te) = discreternd(nobs_te, prob);
  Ytest(inds_te) = randn(nobs_te,1)./phis(Ztest(inds_te))' + means(Ztest(inds_te))';
  Xtest(inds_te) = times(t).*ones(nobs_te,1);
  Xtest_inds(inds_te) = t;

end

Ycell = cell(1,num_times);
Ytest_cell = cell(1,num_times);
for i = 1:num_times
  Ycell{i} = Y(X == times(i))';
  Ytest_cell{i} = Ytest(Xtest == times(i))';
end
% Determine regions for sngp
