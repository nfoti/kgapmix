function [] = list_params()
% LIST_PARAMS Print list of parameters that can passed to params struct for
%  kgap_gmm_finite.

fprintf(' Required arguments:\n');
fprintf('   nsamp : number of samples\n');
fprintf('   S : initial assignments of data to atoms (vector)\n');
fprintf('   Pi : initial atom mass (vector)\n');
fprintf('   theta : initial mixture centers (vector)\n');
fprintf('   phi : initial mixture precisions (vector)\n');
fprintf('   mus : initial dictionary of covariate locations (matrix)\n');
fprintf('   psis : initial dictionary of kernel precisions (vector)\n');
fprintf('   mu_inds : initial covariate location of each atom (vector)\n');
fprintf('   psi_inds : initial kernel precision of each atom (vector)\n');
fprintf('   kernfun : kernel function to use (function handle)\n');

fprintf('\n');

fprintf(' Optional arguments (with defaults denoted in [])\n');
fprintf('   thin [1]: number of samples to skip when collecting\n');
fprintf('   thmean0 [0]: prior mean of thetas (scalar)\n');
fprintf('   thtau0 [1]: prior precision of thetas (scalar)\n');
fprintf('   a0 [1/K]: shape parameter for gamma prior on Pi (K = #atoms)\n');
fprintf('   b0 [1]: rate parameter for gamma prior on Pi\n');
fprintf('   c0 [1]: shape parameter for gamma prior on phis\n');
fprintf('   d0 [1]: rate parameter for gamma prior on phis\n');
fprintf('   U [uniform]: prior on entries of mus (if scalar assumed uniform\n');
fprintf('   V [uniform]: prior on entries of psis (if scalar assumed uniform\n');
fprintf('   sampleMu [false]: whether to sample covariate locations for atoms\n');
fprintf('   samplePhi [false]: whether to sample kernel scale for atoms\n');
fprintf('   sampleVg [false]: whether to sample Vg\n');
fprintf('   sampleVstar [false]: whether to sample Vstar\n');
fprintf('   sampleAlpha [false]: whether to sample alpha\n');
fprintf('   samplePhi [true]: whether to sample cluster precisions\n');
fprintf('   singlePhi [false]: indicates if all components use same precision\n');
fprintf('   ispreds [false] : whether to predict points\n');
fprintf('   predvals : vector of parameters to estimate predictive distributio of\n');
fprintf('   Xpred : matrix of covariates to compute a predictive distribution for\n');

fprintf(' Miscellaneous\n');
fprintf('   rng_type [mt19937ar]: type of random number generator\n');
fprintf('   rng_seed : seed for random stream for sampler (to pick up where we left off)\n');
