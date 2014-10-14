% Run slice sampler on CMB data with square exponential kernel.
%

sdir = regexp(pwd, '/', 'split');
if ~strcmp(sdir{end},'experiments')
  fprintf('Run this script in experiments/ directory\n');
  return;
end

addpath ../gmm/slice;
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
  
  mus = unique(MMtr);
  psis = 1/.005; % precision of 200
  
  % Initialize parameters - These settings from Favaro and Teh
  vbar = var(TTtr);
  thmean0 = mean(TTtr);
  % This is to be consistent with sngp sampler
  thtau0 = .5/vbar;
  % These are determined from the data using the fact that at a covariate value
  % the variance seems to be about 0.02.
  % Might have to play with these a bit
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
  
  samplePhi = true;
  singlePhi = false;
  phi = 1/0.09; % Determined from the data
  
  nburn = 10000;
  nsamp = 5000;
  thin = 1;
  
  ispreds = true;
  [XXun,~,~] = unique(times,'rows');
  npreds = 200;
  predvals = [TTte' linspace(min(TT)-2,max(TT)+2,npreds)];
  Xpred = XXun;
  
  init = init_params_struct(TTtr, XXtr, 'nburn', nburn, 'nsamp', nsamp, 'thin', thin, ...
                                 'mus', mus, 'psis', psis, ...
                                 'thmean0', thmean0, 'thtau0', thtau0, ...
                                 'c0', c0, 'd0', d0, 'L', L, ...
                                 'alpha', alpha, 'phi', phi, ...
                                 'sampleVg', sampleVg, 'sampleVstar', sampleVstar, ...
                                 'sampleAlpha', sampleAlpha, ...
                                 'samplePhi', samplePhi, 'singlePhi', singlePhi, ...
                                 'U', U, 'V', V, 'kernfun', @sekern, ...
                                 'ispreds', ispreds, ...
                                 'Xpred', Xpred, 'predvals', predvals ...
                                 );
  
  params = kgap_gmm_slice(TTtr,XXtr,init);
                    
  % Evaluate
  
  fprintf('Motorcycle se slice:\n');
  
  ittime = params.ittime;
  fprintf('Avg. time / 100 iters: %.2f (+- %.2f)\n', mean(ittime), 1.96*std(ittime)/sqrt(numel(ittime)));
  avg100_slice = mean(ittime);
  
  predinds = 1:ntest;
  pinds = sub2ind(size(params.lpred), te_tind', predinds);
  pll_slice_rb = params.lpred(pinds); % log-lik of held-out points
  fprintf('Avg. RB pred. log-lik: %.2f (+- %.2f)\n', mean(pll_slice_rb), 1.96*std(pll_slice_rb)/sqrt(numel(pll_slice_rb)));
  
  [pll_slice,~,~] = pred_llik_slice(TTte, XXte, params);
  fprintf('Avg. pred. log-lik: %.2f (+- %.2f)\n', mean(pll_slice), 1.96*std(pll_slice)/sqrt(numel(pll_slice)));
  
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
  predvals = linspace(min(TT)-2,max(TT)+2,npreds);
  
  system('mkdir -p motor_se_slice');
  save(['motor_se_slice/motor_se_slice_' num2str(si) '.mat'], 'pll_slice', 'pll_slice_rb', 'pl_slice', 'plse_slice', ...
       'avg100_slice','t_hat_slice', 'predvals', 'XXind', 'KK', 'Ls', 'alphas', 'Ss', 'params');
end

rmpath ../gmm/slice;
rmpath ../util;
