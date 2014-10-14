% Run sngp sampler on CMB data.
%

sdir = regexp(pwd, '/', 'split');
if ~strcmp(sdir{end},'experiments')
  fprintf('Run this script in experiments/ directory\n');
  return;
end

addpath ../sngp_basic;
addpath ../util;

load nips_seeds.mat;

for si = 1:1
%for si = 2:numel(seeds)

  s = RandStream('mt19937ar', 'Seed', seeds(si));
  RandStream.setGlobalStream(s);
  
  % Loads cmb data and takes first 600 data points (they're ordered)
  load_cmb600;
  N = numel(TT);
  pd = pdist2(MM, times);
  [~,tind] = min(pd,[],2);
  
  rp = randperm(N);
  ntest = floor(.1*N);
  ntrain = N - ntest;
  trinds = rp(1:ntrain);
  teinds = rp((ntrain+1):end);
  
  MMtr = MM(trinds);
  TTtr = TT(trinds);
  MMte = MM(teinds);
  TTte = TT(teinds);
  tr_tind = tind(trinds);
  te_tind = tind(teinds);
  
  TTtrain_cell = cell(1,numel(times));
  TTtest_cell = cell(1,numel(times));
  for i = 1:numel(times)
    TTtrain_cell{i} = TTtr(tr_tind==i)';
    TTtest_cell{i} = TTte(te_tind==i)';
  end
  
  % Set parameters
  burnin = 2000;
  num_iter = 5000;
  
  meas = 0.2;
  
  var_x = 0.02; % Determined from the data
  var_means = 1;
  kern_tau = .1;
  
  [rslt, diagn, KK, ittime] = mcmc_silent(TTtrain_cell, times, kern_tau, burnin, ...
                                          num_iter, var_x, var_means, meas);
  
  % Evaluate
  
  fprintf('Avg. time / 100 iters: %.2f (+- %.2f)\n', mean(ittime), 1.96*std(ittime)/sqrt(numel(ittime)));
  avg100_finite = mean(ittime);
  
  [pll_sngp, pl_sngp, pl2_sngp, ~,~] = pred_llik_sngp(TTtest_cell, rslt, var_x);
  fprintf('Avg. pred. log-lik: %.2f (+- %.2f)\n', mean(pll_sngp/num_iter), 1.96*std(pll_sngp/num_iter)/sqrt(numel(pll_sngp)));
  
  % Compute predictive density for a dense sequence of points at all times (for figures)
  predvals = linspace(min(TTte)-2,max(TTte)+2,200);
  TTtest_lsp_cell = cell(1,numel(TTtest_cell));
  for i = 1:numel(TTtest_cell)
    TTtest_lsp_cell{i} = predvals;
  end
  [~,pl_lsp_sngp,pl2_lsp_sngp,~,XXind] = pred_llik_sngp(TTtest_lsp_cell, rslt, var_x);
  
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
  
  % Save output
  system('mkdir -p cmb_box_sngp');
  save(['cmb_box_sngp/cmb_box_sngp_' num2str(si) '.mat'], 'pll_sngp', 'pl_sngp', 'pl2_sngp', 't_hat_sngp', ...
       'avg100_finite', 'predvals', 'XXind','pl_lsp_sngp','pl2_lsp_sngp', 'KK', 'times', 'rslt');
end

rmpath ../sngp_basic;
rmpath ../util;
