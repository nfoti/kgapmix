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

validate_int_pos = @(n)isscalar(n) && n > 0 && n == round(n);
validate_int_nz = @(n)isscalar(n) && n >= 0 && n == round(n);
validate_vec = @(v)isvector(v);
validate_vec_pos = @(v)isvector(v) && all(v >= 0);
validate_matrix = @(m)ismatrix(m);
validate_ind_vec = @(v)isvector(v) && all(v > 0 & v == round(v));
validate_float = @(x)isscalar(x);
validate_float_pos = @(x)isscalar(x) && x > 0;

ip.addParamValue('nburn', 0, validate_int_nz);
ip.addParamValue('nsamp', 1, validate_int_pos);
ip.addParamValue('mus', -1, validate_matrix);
ip.addParamValue('psis', -1, validate_vec_pos);
ip.addParamValue('kernfun', -1);

ip.addParamValue('thin', 1, validate_int_pos);
ip.addParamValue('thmean0', 0, validate_float);
ip.addParamValue('thtau0', 1, validate_float_pos);
ip.addParamValue('c0', 1, validate_float_pos);
ip.addParamValue('d0', 1, validate_float_pos);
ip.addParamValue('L', 10^-2, validate_float_pos);
ip.addParamValue('alpha', 1, validate_float_pos);

ip.addParamValue('phi', 1, validate_float_pos);

ip.addParamValue('U', -1, @(u)(validate_float_pos(u)||validate_vec_pos(u)));
ip.addParamValue('V', -1, @(u)(validate_float_pos(u)||validate_vec_pos(u)));
ip.addParamValue('sampleMu', false, @(b)islogical(b));
ip.addParamValue('samplePsi', false, @(b)islogical(b));
ip.addParamValue('sampleVg', false, @(b)islogical(b));
ip.addParamValue('sampleVstar', false, @(b)islogical(b));
ip.addParamValue('sampleAlpha', false, @(b)islogical(b));
ip.addParamValue('samplePhi', true, @(b)islogical(b));
ip.addParamValue('singlePhi', false, @(b)islogical(b));

ip.addParamValue('ispreds', false, @(b)islogical(b));
ip.addParamValue('predvals', Inf, validate_vec);
ip.addParamValue('Xpred', Inf, validate_matrix);

validrng_types = {'mcg16807','mlfg6331_64','mrg32k3a','mt19937ar', ...
                  'shr3cong','swb2712'};
validate_rng_type = @(x)validatestring(x,validrng_types);
ip.addParamValue('rng_type', 'mt19937ar', validate_rng_type);
ip.addParamValue('rng_seed', -1, validate_int_pos);

ip.parse(varargin{:});
init = ip.Results;

% Validate initial settings
if init.mus == -1
  error(err_id, 'mus is required');
end
if init.psis == -1
  error(err_id, 'psis is required');
end
if isnumeric(init.kernfun)
  error(err_id, 'kernfun is required');
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

if init.ispreds
  if isscalar(init.predvals) && isinf(init.predvals)
    error('Must specify a range of parameter values to predict parameters at');
  end
  if isscalar(init.Xpred) && isinf(init.Xpred)
    error('Must specify a range of covariates to compute predictive distributions for')';
  end
end

