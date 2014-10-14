function [params] = kgap_gmm_slice(Y,X,params)
% KGAP_FINITE Kernel GaP Gaussian mixture model using truncated 
% representation and square exponential kernel.
%
% [] = kgap_finite()
%
% Inputs:
%   Y :  N vector of observations
%   X : [Nxp] covariate matrix (covariates live in p dimensions)
%
% Returns:
%
%
%

if params.rng_seed >= 0
  s = RandStream(params.rng_type,'Seed', params.rng_seed);
  RandStream.setGlobalStream(s);
end

%% Notes
% Note: Below we assume there are K represented atoms
% L : Smallest weighted atom size allowed to be sampled
% S : [Nx1] vector with cluster labels
% theta : K vector with means of mixture components
% phi : K vector of precisions of mixture components
% Pi : K vector of atom weights
% KK : [GxK] matrix with K(xstar_g, pi_k) in (g,k) entry
% V : G-vector with v_g = sum( K(xstar_g,mu)*Pi )
% mus : [Pxp] matrix that is dictionary of unique covariate values
% mu_inds : K vector that indexes into rows of mus dictionary
% psis : Q vector of unique dispersions
% psi_inds : K vector that indexes into psis dictionary
% alpha : mass parameter of base measure on Theta (mass of covariate
%         measure assumed to be 1 for easier sampling)

%% Initialize sampler

thtau0 = params.thtau0;
thmean0 = params.thmean0;
c0 = params.c0;
d0 = params.d0;
mus = params.mus;
psis = params.psis;
alpha = params.alpha;
U = params.U;  % Prior on kernel locations
%V = params.V; % Prior on kernel widths

N = numel(Y);
Ntot = N;

ispreds = params.ispreds; % boolean if we have to predict points
if ispreds == 1
  predvals = params.predvals;       % Values to compute pred. dist. of
  Xpred = params.Xpred;             % Covariate values to compute pred. dist. at
  Npred = size(Xpred,1);
  pred = zeros(Npred, numel(predvals)); % Store the pred. dist.
  pred2 = zeros(Npred, numel(predvals)); % pred. sum of squares for variance
  lpred = zeros(Npred, numel(predvals)); % predictive log-likelihood
  lpred2 = zeros(Npred, numel(predvals)); % predictive log-likelihood squares
  Ntot = Ntot + Npred;
else
  ispreds = 0;
  Xpred = [];
end


if numel(psis) ~= 1
  error('Only can have a single kernel width for now');
end

% find unique covariates and make map from data points to unique covariate
if ispreds
  [xstar,~,y2xstar] = unique([X ; Xpred], 'rows');
else
  [xstar,~,y2xstar] = unique(X,'rows');
end
nn = accumarray(y2xstar, 1); % number of data points assigned to each unique covariate (including predictive covariates)

S = (1:Ntot)';

L = params.L; % 10^-2 reasonable
K = Ntot;

% We may want to only sample Vg for the observed covariates, but for now let's
% just do it for all of them (prediction covariates will just have a sample
% size of 1).
G = size(xstar,1);

% Set up initial atoms
Pi = zeros(1, K);
for i = 1:K
  % Rejection sample from a truncated (from below) gamma distribution
  Pi(i) = draw_DP(L);
end
theta = randn(1,K)./sqrt(thtau0) + thmean0;

if params.samplePhi
  if isfield(params,'phi')
    if numel(params.phi) ~= K
      phi = params.phi.*ones(1,K);
    else
      phi = params.phi;
    end
  elseif params.singlePhi
    phi = gamrnd(c0, 1/d0).*ones(1,K);
  else
    phi = gamrnd(c0, 1/d0, 1, K);
  end
else
  phi = params.phi(1).*ones(1,K);
end

tmp = pdist2([X ; Xpred], mus);
[~,mu_inds] = min(tmp,[],2);
mu_inds = mu_inds';

psi_inds = ones(1,K); % Set atom dispersion to largest value

% Construct kernel values between data points and atom locations
KK = params.kernfun(xstar, mus, mu_inds, psis, psi_inds);

% Set up slice variables and v's
KKmu = bsxfun(@times, KK, Pi);
Vg = gamrnd(nn, 1./(sum(KKmu,2)));
inds = sub2ind([G K], y2xstar, S);
tmp = KKmu(inds);
uu = rand(1,Ntot).*(tmp(:)' - L) + L; % Uniform on [L, (KK.*Pi)(y2xstar,S)] using numpy-type indexing in the latter which is why we need sub2ind above
L = min(uu); % Update L
Vstar = sum(Vg);

nburn = params.nburn;
nsamp = params.nsamp;
thin = params.thin;
niter = nburn + nsamp*thin;

S_samp = zeros(Ntot,nsamp);
Pi_samp = cell(1,nsamp);
Theta_samp = cell(1,nsamp);
Phi_samp = cell(1,nsamp);
mu_inds_samp = cell(1,nsamp);
psi_inds_samp = cell(1,nsamp);
alpha_samp = zeros(1,nsamp);

% Debug collections
Vg_samp = zeros(G,nsamp);
uu_samp = zeros(Ntot,nsamp);
L_samp = zeros(1,nsamp);
Vstar_samp = zeros(1,nsamp);

Vg_accept = zeros(G,1);
Vg_count = zeros(G,1);
Vstar_accept = 0;
Vstar_count = 0;

P = size(mus,1);
%Q = numel(psis);

ittime = zeros(1,nsamp);
itind = 1;

%% Run sampler

sidx = 1;
for iter = 1:niter

  if mod(iter, 100) == 1
    start_t = tic;
  end
  
  %% Sample Vg
  if params.sampleVg
    for g = 1:G
      Jg = KK(g,:).*Pi;
      vgnew = Vg(g) * exp(0.5*randn); % Same as random walk on log-scale, 0.5 is step size

      llnew = (nn(g))*log(vgnew) - vgnew*sum(Jg);  % Check first term, Griffin has nn(g) (probably from jacobian for doing it on log-space)

      % MAKE SURE DIMENSIONS ARE OK
      % Estimate "kernelized tail"
      % IGNORE "kernelized" tail for now as it won't contribute much, we
      % can always upper bound it and overestimate it as well
      %{
      ktail = J < L; % Can use this below as well
      tot = 0;
      for p = 1:size(ktail,2)
        vt = Vg;
        vt(g) = vgnew;
        vt(~ktail(:,p)) = 0;
        tot = tot + exp(-vt*J(:,p));
      end
      tot = tot / nnz(J);
      llnew = llnew + log(tot);
      %}

      % Compute tail of large process (these atoms so small don't contribute
      % anywhere)
      vt = Vg;
      vt(g) = vgnew;
      KKmu = params.kernfun(xstar, mus, 1:P, psis, ones(1,P));
      tmp = vt'*KKmu;
      tot = sum( -log(1+tmp) - expint(L*(1+tmp)) );
      tot = tot / P;
      llnew = llnew + alpha*tot;

      % Do the same for the old value of Vg(g)
      vg = Vg(g);
      llold = (nn(g))*log(vg) - vg*sum(Jg);

      % Estimate "kernelized tail"
      %{
      tot = 0;
      for p = 1:size(ktail,2)
        vt = Vg;
        vt(~ktail(:,p)) = 0;
        tot = tot + exp(-vt*J(:,p));
      end
      tot = tot / nnz(J);
      llold = llold + log(tot);
      %}

      % Compute raw tail
      KKmu = params.kernfun(xstar, mus, 1:P, psis, ones(1,P));
      tmp = Vg'*KKmu;
      tot = sum( -log(1+tmp) - expint(L*(1+tmp)) );
      %for p = 1:P
      %  KKmu = params.kernfun(xstar, mus, p, psis, 1); % Have to do this in case this location isn't used by a Pi
      %  tmp = Vg'*KKmu;
      %  % Check 2nd '-' sign here, Griffin says what's here but seems like it should be '+'
      %  tot = tot - log(1+tmp) - expint(L*(1+tmp));
      %end
      tot = tot / P;
      llold = llold + alpha*tot;

      % MH
      lacc = llnew - llold;
      accept = 1;
      if lacc < 0
        accept = exp(lacc);
      end

      Vg_count(g) = Vg_count(g) + 1;
      % Griffin unconditionally adds accept to V_accept, may be useful rather
      % than below
      if rand < accept
        Vg(g) = vgnew;
        Vg_accept(g) = Vg_accept(g) + 1;
      end

    end
  end
  
  %% Sample Vstar
  % Only sample Vstar if have > 1 unique covariates
  if G > 1 && params.sampleVstar
    B = Vg ./ Vstar;
    J = bsxfun(@times, KK, Pi);
    vsnew = Vstar * exp(0.5*randn); % Same as random walk on log-scale, 0.5 is step size

    llnew = (Ntot)*log(vsnew) - vsnew*sum(B'*sum(J,2));

    % kernelized tail for NEW Vstar
    % ignored for now

    % raw tail for NEW Vstar
    tot = 0;
    KKmu = params.kernfun(xstar, mus, 1:P, psis, ones(1,P));
    tmp = vsnew*B'*KKmu;
    tot = sum( -log(1+tmp) - expint(L*(1+tmp)) );
    tot = tot / P;
    llnew = llnew + alpha*tot;

    llold = (Ntot)*log(Vstar) - Vstar*sum(B'*sum(J,2));

    % kernelized tail for OLD Vstar
    % Ignored for now

    % raw tail for OLD Vstar
    KKmu = params.kernfun(xstar, mus, 1:P, psis, ones(1,P));
    tmp = Vstar*B'*KKmu;
    tot = sum( -log(1+tmp) - expint(L*(1+tmp)) );
    tot = tot / P;
    llold = llold + alpha*tot;

    % MH
    lacc = llnew - llold;
    accept = 1;
    if lacc < 0
      accept = exp(lacc);
    end

    Vstar_count = Vstar_count + 1;
    % Griffin unconditionally adds accept to V_accept, may be useful rather
    % than below
    if rand < accept
      Vstar = vsnew;
      Vstar_accept = Vstar_accept + 1;
    end
  
  end
  
  %% Sample jumps with data allocated (remove jumps with no data right after)
  nk = accumarray(S, 1, [K 1]);
  holds = [];
  counts = zeros(K,1);
  counter = 0;
  for k = 1:K
    if nk(k) > 0
      holds = [holds k]; %cat(2, holds, k);
      counter = counter + 1; % These two lines make the map from old inds to new inds
      counts(k) = counter;
      Pi(k) = gamrnd(nk(k), 1/(1 + Vg'*KK(:,k)));
    end
  end
  
  Pi = Pi(holds);
  KK = KK(:,holds);
  theta = theta(holds);
  phi = phi(holds);
  mu_inds = mu_inds(holds);
  psi_inds = psi_inds(holds);
  S = counts(S); % Map old atom inds to new atom inds
  K = numel(Pi);
  
  %% Sample parameters of remaining used atoms
  % Parameter values for each atom
  Sobs = S(1:N);
  for k = 1:K
  
    kind = Sobs == k;
    nk = sum(kind);
    Yk = Y(kind);
    
    thtau = thtau0 + nk*phi(k);
    thmean = (thtau0*thmean0 + phi(k)*sum(Yk))/thtau;
    
    theta(k) = randn/sqrt(thtau) + thmean;
    
    if params.samplePhi && ~params.singlePhi
      c = c0 + 0.5*nk;
      d = d0 + 0.5*sum( (Yk - theta(k)).^2 );
      phi(k) = gamrnd(c,1./d);
    end
    
  end
  if params.samplePhi && params.singlePhi
    c = c0 + 0.5*N;
    d = d0 + 0.5*sum( (Y - theta(Sobs)').^2 );
    phi = gamrnd(c,1./d).*ones(1,K);
  end
  
  % Update atom covariate locations (we assume constant-width kernels, but
  % put psi_inds(k) below for generality, we could extend the sampler in
  % the future).
  
  % Need this quantity for below
  nngk = zeros(K,G);
  for g = 1:G
    nngk(:,g) = accumarray(S(y2xstar==g), 1, [K 1]);
  end
  
  for k = 1:K
    mu_old = mu_inds(k);
    KKp = params.kernfun(xstar, mus, 1:P, psis, psi_inds(k).*ones(1,P));
    ps = -Pi(k)*(Vg'*KKp) + nngk(k,:)*log(KKp+eps);
    ps = ps - max(ps);
    ps = exp(ps);
    ps = ps ./ sum(ps);
    mu_inds(k) = discreternd(1, ps);
    % Update the value with the new random sample
    if mu_inds(k) ~= mu_old
      KK(:,k) = KKp(:,mu_inds(k));
      %KK(:,k) = params.kernfun(xstar, mus, mu_inds(k), psis, psi_inds(k));
    end
  end
  
  %% Sample slice variables
  KKmu = bsxfun(@times, KK, Pi);
  inds = sub2ind([G K], y2xstar, S);
  tmp = KKmu(inds);
  uu = rand(1,Ntot).*tmp(:)';
  L = min(uu);
  
  %% Sample jump process for unused atoms
  % Mean of Poisson process
  KKmu = params.kernfun(xstar, mus, 1:P, psis, ones(1,P));
  lam = sum( expint(L*(1+Vg'*KKmu)) );
  lam = lam / P;
  
  % Add new atoms and update necessary parameters
  Kadd = min(10^7, poissrnd(alpha*lam)); % 10^7 from Griffin, probably too much now
  Knew = K + Kadd;
  Pi = [Pi zeros(1,Kadd)];
  theta = [theta thmean0+randn(1,Kadd)./sqrt(thtau0)];
  phi = [phi gamrnd(c0, 1/d0, 1, Kadd)];
  mu_inds = [mu_inds randsample(P, Kadd, 1, U)'];
  psi_inds = [psi_inds ones(1,Kadd)];
  KKnew = zeros(G,Knew);
  KKnew(:,1:K) = KK;
  KKnew(:,(K+1):Knew) = params.kernfun(xstar, mus, mu_inds((K+1):end), psis, psi_inds((K+1):end));
  KK = KKnew;
  
  % Sample jump sizes for new atoms
  for k = (K+1):(Knew)
    vv = Vg'*KK(:,k);
    jj = draw_DP( (1+vv)*L ); % This is a convenience thing for the draw_DP function
    Pi(k) = jj / (1+vv);
  end
  K = Knew;
  
  %% Sample cluster indicators
  for i = 1:N
    ps = (uu(i) < KK(y2xstar(i),:).*Pi) .* (exp(-0.5*phi.*(Y(i)-theta).^2) .* sqrt(phi));
    if all(ps < eps)
      keyboard;
    end
    ps = ps ./ sum(ps);
    S(i) = discreternd(1,ps);
  end
  if ispreds
    for i = (N+1):Ntot
      ps = (uu(i) < KK(y2xstar(i),:).*Pi);
      if all(ps < eps)
        keyboard;
      end
      ps = ps ./ sum(ps);
      S(i) = discreternd(1,ps);
    end
  end
  
  %% Sample alpha (concentration parameter)
  if params.sampleAlpha
    aa = 1 + K;
    KKmu = params.kernfun(xstar, mus, 1:P, psis, ones(1,P));
    vv = Vg'*KKmu;
    tmp = sum( log(1+vv) + expint((1+vv)*L) );
    bb = 1 + tmp / P;

    alpha = gamrnd(aa, 1./bb);
  end

  if mod(iter,100) == 0
    ittime(itind) = toc(start_t);
    itind = itind + 1;
  end
  
  %% Grab output
  if iter > nburn && mod(iter,thin) == 0
    S_samp(:,sidx) = S;
    Pi_samp{sidx} = Pi;
    Theta_samp{sidx} = theta;
    Phi_samp{sidx} = phi;
    mu_inds_samp{sidx} = mu_inds;
    psi_inds_samp{sidx} = psi_inds;
    alpha_samp(sidx) = alpha;

    % Rao-Blackwellized estimate of predictive density
    if ispreds
      for pn = 1:Npred
        for j = 1:numel(predvals)
          Jg = KK(y2xstar(N+pn),:).*Pi;
          %prob = (uu(N+pn) < Jg) .* (exp(-0.5*phi.*(predvals(j)-theta).^2) .* sqrt(phi));
          prob = (uu(N+pn) < Jg) .* normpdf(predvals(j), theta, 1./sqrt(phi));
          tmp = sum(prob)/sum(uu(N+pn)<Jg);
          pred(pn,j) = pred(pn,j) + tmp;
          pred2(pn,j) = pred2(pn,j) + tmp.^2;

          ltmp = log(tmp);
          lpred(pn,j) = lpred(pn,j) + ltmp;
          lpred2(pn,j) = lpred2(pn,j) + ltmp.^2;
        end
      end
    end

    % Debug collections
    Vg_samp(:,sidx) = Vg;
    uu_samp(:,sidx) = uu;
    L_samp(sidx) = L;
    Vstar_samp(sidx) = Vstar;
    
    sidx = sidx + 1;
  end
  
  %% Print progress
  if mod(iter,100) == 0
    % Maybe print some sort of log-likelihood
    fprintf('Iter: %d / %d, K: %d (used %d), Kadd: %d, L: %f alpha: %f\n', iter, niter, K, numel(unique(S)), Kadd, L, alpha);
  end
  
end

%% Do accounting to return data
% Pass collected samples out using params. struct
params.S_samp = S_samp;
params.Pi_samp = Pi_samp;
params.Theta_samp = Theta_samp;
params.Phi_samp = Phi_samp;
params.mu_inds_samp = mu_inds_samp;
params.psi_inds_samp = psi_inds_samp;
params.alpha_samp = alpha_samp;

params.Vg_samp = Vg_samp;
params.Vstar_samp = Vstar_samp;
params.uu_samp = uu_samp;
params.L_samp = L_samp;

params.ittime = ittime(1:(itind-1));

params.Vg_accept = Vg_accept;
params.Vg_count = Vg_count;
params.Vstar_accept = Vstar_accept;
params.Vstar_count = Vstar_count;

% Also pass out current values of variables for easy access
params.S = S;
params.Pi = Pi;
params.theta = theta;
params.phi = phi;
params.mu_inds = mu_inds;
params.psi_inds = psi_inds;
params.alpha = alpha;
params.Vg = Vg;
params.Vstar = Vstar;
params.uu = uu;
params.L = L;

% Output predictive distribution
if ispreds
  predbar = pred / nsamp;
  nm1 = nsamp - 1;
  % unbiased estimate of variance using sufficient statistics
  predvar = (nsamp/nm1).*predbar.^2 - 2*predbar.*pred/nm1 + pred2/nm1;
  params.pred = predbar;
  params.predse = sqrt(predvar / niter);

  lpredbar = lpred / niter;
  lpredvar = (nsamp/nm1).*lpredbar.^2 - 2/nm1*lpredbar.*lpred + lpred2/nm1;
  params.lpred = lpredbar;
  params.lpredse = sqrt(lpredvar / niter);
end
