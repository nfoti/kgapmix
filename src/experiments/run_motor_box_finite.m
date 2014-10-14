% Run finite smapler on motorcycle data with box kernel.
%

sdir = regexp(pwd, '/', 'split');
if ~strcmp(sdir{end},'experiments')
  fprintf('Run this script in experiments/ directory\n');
  return;
end

addpath ../gmm/finite;
addpath ../util;

load nips_seeds.mat;

%for si = 1:1
for si = 2:numel(seeds)

  s = RandStream('mt19937ar', 'Seed', seeds(si));
  RandStream.setGlobalStream(s);
  
  % Loads motorcycle data
  load ~/work/data/motor.mat;
  time = time-min(time);
  time = time./max(time);
  % Just renamed the independent and response vars for ease of porting code
  MM = time; 
  TT = (acc-mean(acc))./std(acc);
  clear acc description resvar strata time;
  
  times = 0.1*(.5:.5:9)';
  N = numel(TT);
  pd = pdist2(MM, times);
  [~,tind] = min(pd,[],2);
  
  rp = randperm(N);
  ntest = floor(.1*N);
  ntrain = N - ntest;
  N = ntrain;
  trinds = rp(1:ntrain);
  teinds = rp((ntrain+1):end);
  
  MMtr = MM(trinds);
  TTtr = TT(trinds);
  MMte = MM(teinds);
  TTte = TT(teinds);
  tr_tind = tind(trinds);
  te_tind = tind(teinds);
  
  XXtr = times(tr_tind);
  XXte = times(te_tind);
  
  KK = 100; % Determine from number of clusters used by the slice sampler
  a0 = 1/KK;
  b0 = 1;
  
  mus = unique(MMtr);
  psis = .1;
  
  mu_inds = randsample(size(mus,1),KK,true);
  psi_inds = ones(1,KK);
  
  % Initialize parameters - These settings from Favaro and Teh
  vbar = var(TTtr);
  thmean0 = mean(TTtr);
  % This is to be consistent with sngp sampler
  thtau0 = .5/vbar;
  c0 = 2;
  d0 = vbar/5;
  
  P = numel(mus);
  U = (1/P).*ones(1, P);
  V = 1;
  
  Kern = boxkern(XXtr, mus, mu_inds, psis, psi_inds);
  
  S = zeros(N,1);
  for i = 1:N
    S(i) = randsample(KK, 1, true, Kern(i,:)./sum(Kern(i,:)));
  end
  
  Pi = a0/b0.*ones(1,KK);
  theta = randn(1,KK)./sqrt(thtau0) + thmean0;
  
  samplePsi = false;
  
  samplePhi = true;
  singlePhi = false;
  phi = 1/0.09; % Determined from the data
  
  nburn = 1;
  thin = 5000;
  
  init_burn = init_params_struct(TTtr, XXtr, 'nsamp', nburn, 'thin', thin, ...
                                 'S', S, 'Pi', Pi, 'theta', theta, ...
                                 'phi', phi, 'mus', mus, 'psis', psis, ...
                                 'mu_inds', mu_inds, 'psi_inds', psi_inds, ...
                                 'thmean0', thmean0, 'thtau0', thtau0, ...
                                 'a0', a0, 'b0', b0, 'c0', c0, 'd0', d0, ...
                                 'U', U, 'V', V, 'kernfun', @boxkern, ...
                                 'samplePsi', samplePsi, ...
                                 'samplePhi', samplePhi, 'singlePhi', singlePhi ...
                                );
  fprintf('Burn in...\n');
  params = kgap_gmm_finite(TTtr,XXtr,init_burn);
  
  nsamp = 5000;
  thin = 1;
  
  init = init_params_struct(TTtr, XXtr, 'nsamp', nsamp, 'thin', thin, ...
                            'S', params.S, 'Pi', params.Pi, ...
                            'theta', params.theta, 'phi', params.phi, ...
                            'mus', mus, 'psis', psis, ...
                            'mu_inds', params.mu_inds, ...
                            'psi_inds', params.psi_inds, ...
                            'thmean0', thmean0, 'thtau0', thtau0, ...
                            'a0', a0, 'b0', b0, 'c0', d0, 'd0', d0, ...
                            'U', U, 'V', V, 'kernfun', @boxkern, ...
                            'samplePhi', samplePhi, 'singlePhi', singlePhi ...
                           );
  
  fprintf('Sampling...\n');
  params = kgap_gmm_finite(TTtr,XXtr,init);
  
  % Evaluate
  
  fprintf('Motorcycle box finite:\n');
  
  ittime = params.ittime;
  fprintf('Avg. time / 100 iters: %.2f (+- %.2f)\n', mean(ittime), 1.96*std(ittime)/sqrt(numel(ittime)));
  avg100_finite = mean(ittime);
  
  [pll_finite, pl_finite, pl2_finite] = pred_llik_finite(TTte, XXte, params);
  fprintf('Avg. pred. log-lik: %.2f (+- %.2f)\n', mean(pll_finite/nsamp), 1.96*std(pll_finite/nsamp)/sqrt(numel(pll_finite)));
  
  % Compute predictive density for a dense sequence of points at all times (for figures)
  npreds = 200;
  predvals = repmat(linspace(min(TT)-2,max(TT)+2,npreds)',numel(times),1);
  XXpred = [];
  XXind = [];
  for i = 1:numel(times)
    XXpred = [XXpred ; times(i).*ones(npreds,1)];
    XXind = [XXind ; i.*ones(npreds,1)];
  end
  [~,pl_lsp_finite,pl2_lsp_finite] = pred_llik_finite(predvals, XXpred, params);
  
  predvals = linspace(min(TT)-2,max(TT)+2,npreds);
  
  KK = zeros(1,nsamp);
  for i = 1:nsamp
    KK(i) = numel(unique(params.S_samp(:,i)));
  end
  
  % Integrated autocorrelation on number of clusters like in Kalli paper
  [ac lags] = xcorr(KK,'coeff');
  zeroidx = find(lags==0);
  ac = ac((zeroidx+1):end);
  C = find(abs(ac) < 2/sqrt(nsamp),1,'first');
  if isempty(C)
    C = numel(ac);
  end
  ac = ac(1:(C-1));
  
  t_hat_finite = 0.5 + sum(ac);
  fprintf('tau = %.3f\n', t_hat_finite);
  
  Ss = params.S_samp;
  
  system('mkdir -p motor_box_finite');
  save(['motor_box_finite/motor_box_finite_' num2str(si) '.mat'], 'pll_finite', 'pl_finite', 'pl2_finite', ...
       'avg100_finite', 't_hat_finite', 'predvals', 'XXind', 'pl_lsp_finite', ...
       'pl2_lsp_finite', 'KK', 'Ss', 'params');
end

rmpath ../util;
rmpath ../gmm/finite;
