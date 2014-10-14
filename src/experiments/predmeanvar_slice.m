function [mf mfvar] = predmeanvar_slice(fn)

% All variables we need will be created
load(fn);

ntime = max(XXind);
dx = predvals(2)-predvals(1);

mf = zeros(1,ntime);
mfvar = zeros(1,ntime);

for i = 1:ntime
  mf(i) = sum(predvals.*pl_slice(XXind==i)'.*dx);
  mfvar(i) = sum( (predvals-mf(i)).^2.*pl_slice(XXind==i)'.*dx );
end