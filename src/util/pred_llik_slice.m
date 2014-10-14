function [pll pl pl2] = pred_llik_slice(Y,X,params)
% PRED_LLIK_FINITE Compute the predictive log-likelihood and likelihood of 
% held-out data for the slice sampler.
%
% [pll pl pl2] = pred_llik_slice(Y,X,params)

mus = params.mus;
psis = params.psis;
mu_inds_samp = params.mu_inds_samp;
psi_inds_samp = params.psi_inds_samp;
pi_samp = params.Pi_samp;
theta_samp = params.Theta_samp;
phi_samp = params.Phi_samp;

nsamp = size(theta_samp,2);
npred = numel(Y);
pll = zeros(1,npred);
pl = zeros(1,npred);
pl2 = zeros(1,npred);

for it = 1:nsamp
  for j = 1:npred
    Ke = params.kernfun(X(j), mus, mu_inds_samp{it}, psis, psi_inds_samp{it});
%     prob = bsxfun(@times,...
%                   params.kernfun(X(j), mus, mu_inds_samp(:,it), psis, psi_inds_samp(:,it)),...
%                   pi_samp(:,it)');
    prob = bsxfun(@times, Ke, pi_samp{it});
    psum = sum(prob,2);
    if psum < eps
       warning('no predictive ll for data %d on sample %d', j, it);
       prob = Ke;
    end
    prob = bsxfun(@rdivide, prob, psum);
    tmp = sum( prob .* normpdf(Y(j), theta_samp{it}, sqrt(1./phi_samp{it})) );

    pll(j) = pll(j) + log(tmp);
    pl(j) = pl(j) + tmp;
    pl2(j) = pl2(j) + tmp.^2;
  end
end
