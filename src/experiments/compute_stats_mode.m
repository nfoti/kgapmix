function compute_stats_mode(dirpath, nsamp)
% Compute stats for Table in NIPS paper "Slice sampling kernel-weighted
% normalized random measure mixture models".
%
% dirpath : path to directory where matfiles with results are stored
% nsamp : number of samples each algorithm used (THEY ALL NEED TO USE THE
%         SAME NUMBER)
%
% Note: All mat files name variables consistently with a trailing
% identifier.

fnames = dir(dirpath);
fnames = {fnames(:).name};
fnames = fnames(~strcmp(fnames,'.'));
fnames = fnames(~strcmp(fnames,'..'));

% Extract mode (finite, slice or sngp) from filename of first mat file
tmp = regexp(fnames{1}, '_', 'split');
mode = tmp{3};

nfiles = numel(fnames);

avg100 = zeros(1,nfiles);
mean_pll = zeros(1,nfiles);
std_pll = zeros(1,nfiles);
mean_pll_rb = zeros(1,nfiles);
std_pll_rb = zeros(1,nfiles);
meanKK = zeros(1,nfiles);
stdKK = zeros(1,nfiles);
tau = zeros(1,nfiles);

% Load data
for i =1:nfiles
   S = load(fullfile(dirpath,fnames{i}));
   
   % Common variables
   meanKK(i) = mean(S.KK);
   stdKK(i) = std(S.KK);
   
   if strcmp(mode, 'sngp')
       try
           avg100(i) = S.avg100_sngp;
       catch ME
           avg100(i) = S.avg100_finite; % I misnamed this for some of them
       end
       pll = S.pll_sngp / nsamp;
       mean_pll(i) = mean(pll);
       std_pll(i) = std(pll);
   elseif strcmp(mode, 'slice')
       avg100(i) = S.avg100_slice;
       pll = S.pll_slice / nsamp;
       pll_rb = S.pll_slice_rb;  % Don't divide by nsamp b/c rb estimator already did that
       mean_pll(i) = mean(pll);
       std_pll(i) = std(pll);
       mean_pll_rb(i) = mean(pll_rb);
       std_pll_rb(i) = std(pll_rb);
   elseif strcmp(mode, 'finite')
       avg100(i) = S.avg100_finite;
       pll = S.pll_finite / nsamp;
       mean_pll(i) = mean(pll);
       std_pll(i) = std(pll);
   end
   
   % Compute autocorrelation of number of clusters for this run and store
   [ac,lags] = xcorr(S.KK,'coeff');
   zeroidx = find(lags==0);
   ac = ac((zeroidx+1):end);
   C = find(abs(ac) < 2/sqrt(numel(S.KK)),1,'first');
   if isempty(C)
       C = numel(ac);
   end
   ac = ac(1:(C-1));
   tau(i) = 0.5 + sum(ac);
   
end



% Compute stats from combined data above
fprintf('%s:\n', mode);
fprintf('  avg. time / 100 iters: %.2f\n', mean(avg100));
fprintf('  avg. pred. log-lik: %.2f (+- %.2f)\n', ...
   mean(mean_pll), 1.96*std(mean_pll)/sqrt(nfiles));
if strcmp(mode, 'slice')
fprintf('  RB avg. pred. log-lik: %.4f (+- %.4f)\n', ...
   mean(mean_pll_rb), 1.96*std(mean_pll_rb)/sqrt(nfiles));
end
fprintf('  avg. # clusters: %.2f (+- %.2f)\n', ...
    mean(meanKK), 1.96*std(meanKK)/sqrt(nfiles));
fprintf('  tau: %.3f (+- %.3f)\n', mean(tau), 1.96*std(tau)/sqrt(nfiles));
