means = [0 -1 1];
taus = [1/.1 1/.2 1/.1];
p = [.2 .4 .4];%[.5 .25 .25];
mus = [1 ; 3]; % raw covariate locations
mu_inds = [1 ; 2 ; 2]; % index into raw covariate locations

Z = randsample(3, N, true, p);
Z = sort(Z);
Y = randn(1,N)./sqrt(taus(Z)) + means(Z);

% data point location by dereferencing raw location that corresponding mean
% uses
XX = mus(mu_inds(Z));

clear mus mu_inds

% Plot true density
%{
xx = linspace(min(Y)-.5, max(Y)+.5, 100);
dens = zeros(1, numel(xx));
for j = 1:3
  dens = dens + p(j).*normpdf(xx, means(j), 1/sqrt(taus(j)));
end
%}
