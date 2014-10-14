function [params] = kgap_gmm_finite(Y,X,params,name)
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

% Initialize sampler

if params.rng_seed >= 0
  s = RandStream(params.rng_type,'Seed', params.rng_seed);
  RandStream.setGlobalStream(s);
end

% S : [Nx1] vector with cluster labels
% theta : K vector with means of mixture components
% phi : K vector of precisions of mixture components
% Pi : K vector of atom weights
% mus : [Pxp] matrix that is dictionary of unique covariate values
% mu_inds : K vector that indexes into rows of mus dictionary
% psis : K vector of unique dispersions
% psi_inds : K vector that indexes into psis dictionary

N = numel(Y);
K = numel(params.theta);
nsamp = params.nsamp;
thin = params.thin;
niter = nsamp*thin;

% Construct kernel values between data points and atom locations
KK = params.kernfun(X, params.mus, params.mu_inds, params.psis, params.psi_inds);

S_samp = zeros(N,nsamp);
Pi_samp = zeros(K,nsamp);
Theta_samp = zeros(K,nsamp);
Phi_samp = zeros(K,nsamp);
Phi_samp = zeros(K,nsamp);
mu_inds_samp = zeros(K,nsamp);
psi_inds_samp = zeros(K,nsamp);

ittime = zeros(1,nsamp);
itind = 1;

mus = params.mus;
psis = params.psis;
P = size(mus,1);
Q = numel(psis);

% Run sampler

lp = zeros(1,niter);

logf = fopen(['./kgam_gmm_finite' name '.log'], 'w');
fprintf(logf, 'Sampling\n');

sidx = 1;
for iter = 1:niter

  if mod(iter, 10) == 1
    tstart = tic;
  end
  
  % Sample s_i
  lPi = log(params.Pi+eps);
  for i = 1:N
    ps = -Inf.*ones(1,K);
    for k = 1:K
      if KK(i,k) > 0 && params.Pi(k) > 0
        ps(k) = log(KK(i,k)+eps) + lPi(k) + 0.5*log(params.phi(k)) ...
              - 0.5*params.phi(k)*(Y(i)-params.theta(k)).^2;
      end
    end
    ps = ps - max(ps); % For numerical stability
    ps = exp(ps);
    ps = ps ./ sum(ps);
    params.S(i) = discreternd(1,ps);
  end

  % Sample pi_k
  for k = 1:K
    
    inds = sub2ind(size(KK), (1:N)', params.S);
    KKpi = bsxfun(@times, KK, params.Pi);
    denom = sum(KKpi,2);
    numer = KKpi(inds);
    
    uu = rand(N,1).*(numer./denom);
    
    % Find lower truncation
    ki = params.S==k;
    lb = max( uu(ki) ./ (1-uu(ki)+eps) ...
              .* (sum(KKpi(ki,[1:(k-1) (k+1):end]),2) ./ (KK(ki,k)+eps)) );
    if isempty(lb)
      lb = 10^-9;
    end
            
    % Find upper truncation
    nki = params.S~=k;
    ub = min(1./(KK(nki,k)+eps) .* ( KKpi(inds(nki))./(uu(nki)+eps) - sum(KKpi(nki,[1:(k-1) (k+1):end]),2) ));
    if isempty(ub)
      ub = 10^9;
    end
    
    params.Pi(k) = max(gamrndT(1, params.a0, params.b0, lb, ub), 10^-9);
    
  end
  
  % Only sample mu_inds and psi_inds after 
  % Sample mu_k
  if params.sampleMu
    inds = sub2ind(size(KK), (1:N)', params.S);
    for k = 1:K
      mu_old = params.mu_inds(k);
      ps = zeros(1,P);
      KKp = KK;
      for p = 1:P
        KKp(:,k) = params.kernfun(X, mus, p, psis, params.psi_inds(k));
        %KKpi_p = bsxfun(@times, KKp, params.Pi);
        %KKpi_p = bsxfun(@rdivide, KKpi_p, sum(KKpi_p,2));
        ps(p) = log(params.U+eps) + sum( log(KKp(inds)+eps) + log(params.Pi(params.S)) );
        ps(p) = ps(p) - sum( sum(KKp,2)  );
        %ps(p) = log(params.U+eps) + sum( KKpi_p + eps);
      end
      ps = ps - max(ps);
      ps = exp(ps);
      ps = ps ./ sum(ps);
      params.mu_inds(k) = discreternd(1, ps);
      % Update the value with the new random sample
      if params.mu_inds(k) ~= mu_old
        KK(:,k) = params.kernfun(X, mus, params.mu_inds(k), psis, psi_inds(k));
      end
    end
  end
  
  % Sample psi_k
  if params.samplePsi
    for k = 1:K
      psi_old = params.psi_inds(k);
      ps = zeros(1,P);
      KKq = KK;
      for q = 1:Q
        KKq(:,k) = params.kernfun(X, mus, params.mu_inds(k), psis, q);
        KKpi_q = bsxfun(@times, KKq, params.Pi);
        KKpi_q = bsxfun(@rdivide, KKpi_q, sum(KKpi_q,2));
        ps(q) = log(params.V+eps) + sum( KKpi_q + eps);
      end
      ps = ps - max(ps);
      ps = exp(ps);
      ps = ps ./ sum(ps);
      params.psi_inds(k) = discreternd(1, ps);
      % Update the value with the sampled value
      if params.psi_inds(k) ~= psi_old
        KK(:,k) = params.kernfun(X, mus, params.mu_inds(k), psis, psi_inds(k));
      end
    end
  end
  
  % Sample theta_k and phi_k
    for k = 1:K

      kind = params.S == k;
      nn = sum(kind);
      Yk = Y(kind);
      
      thtau0 = params.thtau0;
      thmean0 = params.thmean0;
      thtau = thtau0 + nn*params.phi(k);
      thmean = (thtau0*thmean0 + params.phi(k)*sum(Yk))/thtau;
      
      params.theta(k) = randn*sqrt(1./thtau) + thmean;
      
      if params.samplePhi && ~params.singlePhi
        c = params.c0 + 0.5*nn;
        d = params.d0 + 0.5*sum( (Yk - params.theta(k)).^2 );
        params.phi(k) = gamrnd(c,1./d);
      end
    end
    if params.samplePhi && params.singlePhi
      c = params.c0 + 0.5*N;
      d = params.d0 + 0.5*sum( (Y - params.theta(params.S)').^2 );
      params.phi = gamrnd(c, 1./d).*ones(1,K);
    end

  if mod(iter,10) == 0
    ittime(itind) = toc(tstart);
    itind = itind + 1;
    ittime(itind-1)
  end
  
  lp(iter) = logP(Y,KK,params);
  
  if mod(iter,thin) == 0
    S_samp(:,sidx) = params.S;
    Pi_samp(:,sidx) = params.Pi;
    Theta_samp(:,sidx) = params.theta;
    Phi_samp(:,sidx) = params.phi;
    Phi_samp(:,sidx) = params.phi;
    mu_inds_samp(:,sidx) = params.mu_inds;
    psi_inds_samp(:,sidx) = params.psi_inds;
    
    sidx = sidx + 1;
  end
  
  if mod(iter,1) == 0
    fprintf(logf,'Iter: %d / %d,\n', iter, niter);
  end
  
end

params.S_samp = S_samp;
params.Pi_samp = Pi_samp;
params.Theta_samp = Theta_samp;
params.Phi_samp = Phi_samp;
params.Phi_samp = Phi_samp;
params.mu_inds_samp = mu_inds_samp;
params.psi_inds_samp = psi_inds_samp;
params.lp = lp;

params.ittime = ittime(1:(itind-1));

gs = RandStream.getGlobalStream;
params.rng_seed = gs.Seed;
params.rng_type = gs.Type;

fclose(logf);
