% Generate synthetic data (train and test) and run slice sampler for a box
% kernel.

addpath ../util;
addpath ../gmm/slice;

load nips_seeds.mat;
% This will be a loop

for si = 1:1
%for si = 2:numel(seeds)

  s = RandStream('mt19937ar', 'Seed', seeds(si));
  RandStream.setGlobalStream(s);
  
  % If want to change kernel width do it in function below (it's called kern_tau
  % there)
  gen_fake_data_box;
  
  % This part is sampler specific
  
  mus = unique(X);
  psis = kern_tau;
  
  % Initialize parameters - These settings from Favaro and Teh
  vbar = var(Y);
  thmean0 = mean(Y);
  % This is to be consistent with sngp sampler
  thtau0 = 1./var_means; %.5/vbar;
  c0 = 2;
  d0 = vbar/5;
  
  P = numel(mus);
  U = (1/P).*ones(1, P);
  V = 1;
  
  L = 10^-2;
  
  alpha = 1;
  
  sampleVg = true;
  sampleVstar = true;
  sampleAlpha = true;
  
  samplePhi = false;
  singlePhi = true;
  phi = phis(1);
  
  nburn = 10000;
  nsamp = 5000;
  thin = 1;
  
  %ispreds = false;
  % This is very slow to do
  % Build up one huge vector of points to predict.  Contains test points and
  % complete predictive densities.
  
  ispreds = true;
  [XXun,~,XXi] = unique(X,'rows');
  npreds = 200;
  ntest = numel(Ytest);
  predvals = [Ytest' linspace(min(Y)-5,max(Y)+5,npreds)];
  Xpred = XXun;
  
  
  init = init_params_struct(Y, X, 'nburn', nburn, 'nsamp', nsamp, 'thin', thin, ...
                                 'mus', mus, 'psis', psis, ...
                                 'thmean0', thmean0, 'thtau0', thtau0, ...
                                 'c0', c0, 'd0', d0, 'L', L, ...
                                 'alpha', alpha, 'phi', phi, ...
                                 'sampleVg', sampleVg, 'sampleVstar', sampleVstar, ...
                                 'sampleAlpha', sampleAlpha, ...
                                 'samplePhi', samplePhi, 'singlePhi', singlePhi, ...
                                 'U', U, 'V', V, 'kernfun', @boxkern, ...
                                 'ispreds', ispreds, ...
                                 'Xpred', Xpred, 'predvals', predvals ...
                                 );
                               
  fprintf('Sampling...\n');
  start_t = tic;
  
  params = kgap_gmm_slice(Y,X,init);
  
  elapsed_t = toc(start_t);
  fprintf('Total time: %.2f\n', elapsed_t);
  
  ittime = params.ittime;
  fprintf('Avg. time / 100 iters: %.2f (+- %.2f)\n', mean(ittime), 1.96*std(ittime)/sqrt(numel(ittime)));
  avg100_slice = mean(ittime);
  
  % Evaluate
  predinds = 1:ntest;
  pinds = sub2ind(size(params.lpred), Xtest_inds, predinds);
  pll_slice_rb = params.lpred(pinds); % log-lik of held-out points
  fprintf('Avg. RB pred. log-lik: %.2f (+- %.2f)\n', mean(pll_slice_rb), 1.96*std(pll_slice_rb)/sqrt(numel(pll_slice_rb)));
  
  % Compute the pred. log-like using the corresponding estimator as the sngp and
  % finite
  [pll_slice,~,~,] = pred_llik_slice(Ytest, Xtest, params);
  fprintf('Avg. pred. log-lik: %.2f (+- %.2f)\n', mean(pll_slice/nsamp), 1.96*std(pll_slice/nsamp)/sqrt(numel(pll_slice)));
  
  
  % pred. dens for each covariate (format consistent with sngp version)
  pl_slice = params.pred(:,(ntest+1):end)';
  pl_slice = pl_slice(:);
  plse_slice = params.predse(:,(ntest+1):end)';
  plse_slice = plse_slice(:);
  XXind = [];
  for i = 1:numel(Xpred)
    XXind = [XXind i*ones(1,npreds)];
  end
  
  % Integrated autocorrelation on number of clusters like in Kalli
  KK = zeros(1, nsamp);
  for i = 1:nsamp
    KK(i) = numel(unique(params.S_samp(:,i)));
  end
  [ac lags] = xcorr(KK,'coeff');
  zeroidx = find(lags==0);
  ac = ac((zeroidx+1):end);
  C = find(abs(ac) < 2/sqrt(nsamp),1,'first');
  if isempty(C)
    C = numel(ac);
  end
  ac = ac(1:(C-1));
  
  t_hat_slice = 0.5 + sum(ac);
  fprintf('tau = %.3f\n', t_hat_slice);
  
  Ls = params.L_samp;
  alphas = params.alpha_samp;
  Ss = params.S_samp;
  predvals = linspace(min(Y)-5,max(Y)+5,npreds);
  
  system('mkdir -p synth_box_slice');
  save(['synth_box_slice/synth_box_slice_' num2str(si) '.mat'], 'pll_slice', 'pll_slice_rb', 'pl_slice', ...
       'plse_slice', 'avg100_slice', 't_hat_slice', 'predvals', 'XXind', ...
       'KK', 'Ls', 'alphas', 'Ss', 'params');
end

rmpath ../util;
rmpath ../gmm/slice;
