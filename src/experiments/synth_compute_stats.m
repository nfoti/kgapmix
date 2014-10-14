% Load saved results from three synth_box experiments, compute statistics and
% write data for plots.
%
% TODO: Update to load all 5 result files for each data set

%% sngp
load synth_box_sngp;

fprintf('sngp: Avg. time / 100 iters: %.2f\n', avg100_sngp);

fprintf('sngp: Avg. pred. log-lik: %.2f (+- %.2f)\n', mean(pll_sngp), 1.96*std(pll_sngp)/sqrt(numel(pll_sngp)));

fprintf('sngp: Avg. # clusts: %.2f (+- %.2f)\n', mean(KK), 1.96*std(KK)/numel(KK));

% Integrated autocorrelation on number of clusters like in Kalli paper
[ac,lags] = xcorr(KK,'coeff');
zeroidx = find(lags==0);
ac = ac((zeroidx+1):end);
C = find(abs(ac) < 2/sqrt(numel(KK)),1,'first');
if isempty(C)
  C = numel(ac);
end
ac = ac(1:(C-1));

t_hat_sngp = 0.5 + sum(ac);
fprintf('sngp: tau = %.3f\n', t_hat_sngp);

% Plot predictive density...

clear;


%% slice
load synth_box_slice;

fprintf('slice: Avg. time / 100 iters: %.2f\n', avg100_slice);

fprintf('slice: Avg. pred. log-lik: %.2f (+- %.2f)\n', mean(pll_slice), 1.96*std(pll_slice)/sqrt(numel(pll_slice)));

fprintf('slice: Avg. # clusts: %.2f (+- %.2f)\n', mean(KK), 1.96*std(KK)/numel(KK));

[ac lags] = xcorr(KK,'coeff');
zeroidx = find(lags==0);
ac = ac((zeroidx+1):end);
C = find(abs(ac) < 2/sqrt(numel(KK)),1,'first');
if isempty(C)
  C = numel(ac);
end
ac = ac(1:(C-1));

t_hat_slice = 0.5 + sum(ac);
fprintf('slice: tau = %.3f\n', t_hat_slice);

% Plot predictive density...

clear


%% finite
load synth_box_finite;

fprintf('finite: Avg. time / 100 iters: %.2f\n', avg100_finite);

fprintf('finite: Avg. pred. log-lik: %.2f (+- %.2f)\n', mean(pll_finite), 1.96*std(pll_finite)/sqrt(numel(pll_finite)));

fprintf('finite: Avg. # clusts: %.2f (+- %.2f)\n', mean(KK), 1.96*std(KK)/numel(KK));

[ac lags] = xcorr(KK,'coeff');
zeroidx = find(lags==0);
ac = ac((zeroidx+1):end);
C = find(abs(ac) < 2/sqrt(numel(KK)),1,'first');
if isempty(C)
  C = numel(ac);
end
ac = ac(1:(C-1));

t_hat_finite = 0.5 + sum(ac);
fprintf('finite: tau = %.3f\n', t_hat_finite);

% Plot predictive density...

clear
