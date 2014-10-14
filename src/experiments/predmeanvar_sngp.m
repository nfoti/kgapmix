function [mf mfvar] = predmeanvar_sngp(fn)

% All variables we need will be created
load(fn);

nsamp = numel(KK);

ntime = max(XXind);
dx = predvals(2)-predvals(1);

mf = zeros(1,ntime);
mfvar = zeros(1,ntime);

pl_lsp_sngp = pl_lsp_sngp / nsamp; %#ok

for i = 1:ntime
  mf(i) = sum(predvals.*pl_lsp_sngp(XXind==i).*dx);
  mfvar(i) = sum( (predvals-mf(i)).^2.*pl_lsp_sngp(XXind==i).*dx );
end