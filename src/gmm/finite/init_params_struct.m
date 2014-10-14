function [init] = init_params_struct(Y,X,varargin)
% INIT_PARAMS_STRUCT Set initial values of init structure.  Throws errors
%  for missing fields that are required and does some validation.

err_id = 'init_params_struct:error';

% Before we do anything validate observations
if numel(Y) ~= size(X,1)
  error(err_id,'Differing number of observations and covariates');
end

ip = inputParser;
ip.CaseSensitive = false;
ip.FunctionName = 'init_params_struct';
ip.KeepUnmatched = false;

validate_int_pos = @(n)isnumeric(n) && n > 0 && n == round(n);
validate_vec = @(v)isvector(v);
validate_vec_pos = @(v)isvector(v) && all(v >= 0);
validate_matrix = @(m)ismatrix(m);
validate_ind_vec = @(v)isvector(v) && all(v > 0 & v == round(v));
validate_float = @(x)isscalar(x);
validate_float_pos = @(x)isscalar(x) && x > 0;

ip.addParamValue('nsamp', 1, validate_int_pos);
ip.addParamValue('S', -1, validate_vec);
ip.addParamValue('Pi', -1, validate_vec_pos);
ip.addParamValue('theta', -1, validate_vec);
ip.addParamValue('phi', -1, validate_vec_pos);
ip.addParamValue('mus', -1, validate_matrix);
ip.addParamValue('psis', -1, validate_vec_pos);
ip.addParamValue('mu_inds', -1, validate_ind_vec);
ip.addParamValue('psi_inds', -1, validate_ind_vec);
ip.addParamValue('kernfun', -1);

ip.addParamValue('thin', 1, validate_int_pos);
ip.addParamValue('thmean0', 0, validate_float);
ip.addParamValue('thtau0', 1, validate_float_pos);
ip.addParamValue('a0', -1, validate_float_pos);
ip.addParamValue('b0', 1, validate_float_pos);
ip.addParamValue('c0', 1, validate_float_pos);
ip.addParamValue('d0', 1, validate_float_pos);
ip.addParamValue('U', -1, @(u)(validate_float_pos(u)||validate_vec_pos(u)));
ip.addParamValue('V', -1, @(u)(validate_float_pos(u)||validate_vec_pos(u)));
ip.addParamValue('sampleMu', false, @(b)islogical(b));
ip.addParamValue('samplePsi', false, @(b)islogical(b));
ip.addParamValue('samplePhi', true, @(b)islogical(b));

ip.addParamValue('singlePhi', false, @(b)islogical(b));

validrng_types = {'mcg16807','mlfg6331_64','mrg32k3a','mt19937ar', ...
                  'shr3cong','swb2712'};
validate_rng_type = @(x)validatestring(x,validrng_types);
ip.addParamValue('rng_type', 'mt19937ar', validate_rng_type);
ip.addParamValue('rng_seed', -1, validate_int_pos);

ip.parse(varargin{:});
init = ip.Results;

% Validate initial settings
if init.S == -1
  error(err_id, 'S is required');
end
if init.Pi == -1
  error(err_id, 'Pi is required');
end
if init.theta == -1
  error(err_id, 'theta is required');
end
if init.phi == -1
  error(err_id, 'phi is required');
end
if init.mus == -1
  error(err_id, 'mus is required');
end
if init.psis == -1
  error(err_id, 'psis is required');
end
if init.mu_inds == -1
  error(err_id, 'mu_inds is required');
end
if init.psi_inds == -1
  error(err_id, 'psi_inds is required');
end
if isnumeric(init.kernfun)
  error(err_id, 'kernfun is required');
end

if numel(init.phi) > 1
  szvec = [numel(init.theta) numel(init.phi) numel(init.Pi) ...
           numel(init.mu_inds) numel(init.psi_inds)];
else
  szvec = [numel(init.theta) numel(init.Pi) ...
           numel(init.mu_inds) numel(init.psi_inds)];
end
if numel(unique(szvec)) ~= 1
  fprintf('Number of specified clusters are not consistent:\n');
  fprintf('\t theta: %d\n', numel(init.theta));
  fprintf('\t phi: %d\n', numel(init.phi));
  fprintf('\t Pi: %d\n', numel(init.Pi));
  fprintf('\t mu_inds: %d\n', numel(init.mu_inds));
  fprintf('\t psi_inds: %d\n', numel(init.psi_inds));
  error(err_id, 'Inconsistant number of atoms');
end

K = numel(init.theta);

if init.a0 == -1
  init.a0 = 1/K;
end

if isscalar(init.U)
  init.U = 1/size(init.mus,1);
else
  init.U = bsxfun(@rdivide, init.U, sum(init.U));
end
if isscalar(init.V)
  init.V = 1/numel(init.psis);
else
  init.V = bsxfun(@rdivide, init.V, sum(init.V));
end

if init.thin == 0
  init.thin = 1;
end

% If specified single phi, make sure they're all equal
if init.singlePhi || numel(init.phi == 1)
  init.phi = init.phi(1).*ones(1,K);
end
