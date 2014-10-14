% Run slice sampler on toy time-varying mixture of 3 components

%s = RandStream('mt19937ar', 'Seed', 8675309);
s = RandStream('mt19937ar', 'Seed', floor(rand*1000000));
RandStream.setGlobalStream(s);
fprintf('RandStream Seed: %d\n', RandStream.getGlobalStream.Seed);

% Generate data
N = 500;
Y = [];
XX = [];
Z = [];
means = [];
taus = [];
p = [];
gen_mix3_tv;

spread = max(Y) - min(Y);
mY = mean(Y);
vbar = var(Y);

KK = 10;
a0 = 1/KK;
b0 = 1;

% thmean0 = data mean
% thtau0  = .5/sample variance
% c0 and d0 are chosen so that mean of phi is 10/vbar and variance is
%  50/vbar (taken from Yee Whye's MCMC for NRM paper)
thmean0 = mY;
thtau0 = .5/vbar;
c0 = 2;
d0 = vbar/5;
U = [.5 .5];
V = 1;

L = 10^-2;

alpha = 1;

sampleVg = true;
sampleVstar = true;
sampleAlpha = true;

samplePhi = true;
singlePhi = false;

mus = unique(XX);
psis = 1;

nburn = 10000;
nsamp = 1000;
thin = 5;

ispreds = true;
Xpred = [unique(XX,'rows') ; 0.5 ; 2];
predvals = linspace(min(Y)-0.5, max(Y)+0.5, 50);

init = init_params_struct(Y, XX, 'nburn', nburn, 'nsamp', nsamp, 'thin', thin, ...
                               'mus', mus, 'psis', psis, ...
                               'thmean0', thmean0, 'thtau0', thtau0, ...
                               'c0', c0, 'd0', d0, 'L', L, ...
                               'alpha', alpha, ...
                               'sampleVg', sampleVg, 'sampleVstar', sampleVstar, ...
                               'sampleAlpha', sampleAlpha, ...
                               'singlePhi', singlePhi, 'samplePhi', samplePhi, ...
                               'U', U, 'V', V, 'kernfun', @sekern, ...
                               'ispreds', ispreds, 'Xpred', Xpred, ...
                               'predvals', predvals ...
                              );
fprintf('Sampling...\n');
params = kgap_gmm_slice(Y,XX,init);
