% Generate synthetic data (train and test) and run sngp sampler for a box
% kernel.

sdir = regexp(pwd, '/', 'split');
if ~strcmp(sdir{end},'experiments')
  fprintf('Run this script in experiments/ directory\n');
  return;
end

addpath ../sngp_basic;
addpath ../util;

load nips_seeds.mat;
% We will make this a loop

%for si = 1:1
for si = 1:numel(seeds)

  s = RandStream('mt19937ar', 'Seed', seeds(si));
  RandStream.setGlobalStream(s);
  
  % If want to change kernel width do it in function below (it's called kern_tau
  % there)
  % Creates var_means, kern_tau
  gen_fake_data_box;
  
  % To be comparable with slice sampler
  burnin = 2000;
  num_iter = 5000;
  
  meas = 0.2;
  
  start_t = tic;
  
  [rslt, diagn, KK, ittime] = mcmc_silent(Ycell, times, kern_tau, burnin, ...
                                          num_iter, 1/phis(1), var_means, meas);
  
  elapsed_t = toc(start_t);
  fprintf('Total time: %.2f\n', elapsed_t);
  
  fprintf('Avg. time / 100 iters: %.2f (+- %.2f)\n', mean(ittime), 1.96*std(ittime)/sqrt(numel(ittime)));
  avg100_sngp = mean(ittime);
  
  % Compute predictive likelihood for held out points
  [pll_sngp, pl_sngp, pl2_sngp, ~, ~] = pred_llik_sngp(Ytest_cell, rslt, 1/phis(1));
  fprintf('Avg. pred. log-lik: %.2f (+- %.2f)\n', mean(pll_sngp/nsamp), 1.96*std(pll_sngp/nsamp)/sqrt(numel(pll_sngp)));
  
  % Compute predictive density for a dense sequence of points at all times (for figures)
  predvals = linspace(min(Y)-5,max(Y)+5,200);
  Ytest_lsp_cell = cell(1,numel(Ytest_cell));
  for i = 1:numel(Ytest_cell)
    Ytest_lsp_cell{i} = predvals;
  end
  [~,pl_lsp_sngp,pl2_lsp_sngp,YY,XXind] = pred_llik_sngp(Ytest_lsp_cell, rslt, 1/phis(1));
  %XX = times(XXind);
  
  % Integrated autocorrelation on number of clusters like in Kalli paper
  [ac lags] = xcorr(KK,'coeff');
  zeroidx = find(lags==0);
  ac = ac((zeroidx+1):end);
  C = find(abs(ac) < 2/sqrt(num_iter),1,'first');
  if isempty(C)
    C = numel(ac);
  end
  ac = ac(1:(C-1));
  
  t_hat_sngp = 0.5 + sum(ac);
  fprintf('tau = %.3f\n', t_hat_sngp);
  
  system('mkdir -p synth_box_sngp');
  save(['synth_box_sngp/synth_box_sngp_' num2str(si) '.mat'], 'pll_sngp', 'pl_sngp', 'pl2_sngp', 'avg100_sngp', 't_hat_sngp', ...
       'predvals', 'XXind','pl_lsp_sngp','pl2_lsp_sngp', 'KK', 'rslt');
end

rmpath ../sngp_basic;
rmpath ../util;
