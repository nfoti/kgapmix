function [KK] = boxkern(X, mus, mu_inds, psis, psi_inds)
% BOXKERN Compute box kernel between covariate and atom locations.
%
% [KK] = boxkern(X, mus, mu_inds, psis, psi_inds, PD)
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
%

N = size(X,1);
K = numel(mu_inds);

psi_inds = psi_inds(:)';

KK = bsxfun(@le, pdist2(X,mus(mu_inds,:)), psis(psi_inds));
