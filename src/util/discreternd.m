function r = discreternd(n,p)
% DISCRETERND Sample n random values from discrete distribution p returning 
%   index of sampled value.
%   
%   n : number of random values
%   p : discrete distribution to draw from (should be normalized)
%
%   This code is based on similar code by Finale Doshe-Velez.

if n > 0
  r = ones(n,1);
  % Don't do too much work if only 2 values
  if numel(p) == 2
    for i = 1:n
      if rand < p(2)
        r(i) = 2;
      end
    end
  else
    cump = cumsum(p);
    cump(end) = 1+eps;
    for i = 1:n
      % find first index where cdf is > the random value
      r(i) = find((rand*(1-eps)) < cump, 1); 
    end
  end
else
  r = [];
end
