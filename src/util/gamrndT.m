function g = gamrndT(N,a,b,ll,rr)
% GAMRNDT Sample from truncated gamma random variable.
%
% [g] = gamrndT(N,a,b,ll,rr)
%   N : Number of points to draw
%   a : scale parameter
%   b : rate parameter
%   ll : lower truncation point
%   rr : upper truncation point
%
% Algorithm from "Simulation of truncated gamma variables" by Y. Chung.
% It's basically just doing rejection sampling, but uses a good acceptance
% function depending on the params.

del = rr-ll;
g = zeros(N,1);

x0 = (a-1)/b; % mode of gamma distribution

for i = 1:N

  count = 1;
  while 1 
    x = rand*del + ll;
    num = x^(a-1)*exp(-b*x);
    if x0 >= ll && x0 <= rr && a >= 1
      gx = num / (((a-1)/b)^(a-1)*exp(-(a-1)));
    elseif x0 > rr && a >= 1
      gx = num / (rr^(a-1)*exp(-b*rr));
    elseif (x0 < ll && a >= 1) || a < 1
      gx = num / (ll^(a-1)*exp(-b*ll));
    end
    
    % If we're sampling from something close to a gamma process and count
    % is large, start shrinking the upper bound
    if a < 1 && count > 10
      rr = rr / 2;
      break;
    end
    
    u = rand;
    if u <= gx
      g(i) = x;
      break;
    end
  
    count = count + 1;
    
  end
  
end
