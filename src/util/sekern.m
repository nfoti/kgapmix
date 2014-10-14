function [KK] = sekern(X, mus, mu_inds, psis, psi_inds)
% SEKERN Compute square exponential kernel between covariate locations and atom
%   locations.
%
% [KK] = sekern(X, mus, mu_inds, psis, psi_inds)
%
% Inputs:
%   X : [Nxp] covariate matrix (data in rows)
%   mus : [Pxp] matrix of covariate dictionary
%   mu_inds : K vector indexing into mus dictionary
%   psis : K vector of unique kernel dispersions
%   psi_inds : K vector indexing into psis
%
% Returns:
%   KK : [NxK] matrix of kernel evaluated from x_i to mu_k using psi_k 
%
% Notes:

N = size(X,1);
K = numel(mu_inds);

psi_inds = psi_inds(:)';

KK = exp(-bsxfun(@times, pdist2(X,mus(mu_inds,:)).^2, psis(psi_inds)));
