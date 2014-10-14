function [time2sreg_tab, sreg2time_tab] = build_regions_sngp(t, L)
%
%
% t is a vector of time stamps (should be sorted)
% L is constant window length

% Number of time steps
T = numel(t);
t = t(:)'; % This actually seems faster than reshape on some simple trials

%% Region computations, 
%  This should be changed for other model

% Precompute disjoint regions R that make up the time intervals
tints = zeros(T,2); % First figure out the intervals
for tidx = 1:T
  tints(tidx,:) = [t(tidx) - L , t(tidx) + L];
end

% Compute the regions R
% endlabels indicates if a point was a left end point or a right end point,
% 0 means left, 1 means right.
LEFT = 0;
RIGHT = 1;
R = [];
endpoints = [tints(:,1) ; tints(:,2)];
endlabels = [zeros(size(tints,1),1) ; ones(size(tints,1),1)];
[endpoints perm] = sort(endpoints,'ascend');
endlabels = endlabels(perm);

for tidx = 1:numel(endpoints)-1
  e1 = endpoints(tidx);
  e2 = endpoints(tidx+1);
  t1 = endlabels(tidx);
  t2 = endlabels(tidx+1);
  if t1 == LEFT
    R = [R ; [e1,e2]]; %#ok
  elseif t1 == RIGHT
    if t2 == RIGHT
      R = [R ; [e1,e2]]; %#ok
    else
      ptl = endpoints(endpoints < e1);
      ptlt = endlabels(endpoints < e1);
      if any((e1 - ptl(ptlt == LEFT)) <= L) % Add this interval b/c it sits in another
        R = [R ; [e1,e2]]; %#ok
      end
      % Otherwise we're at a gap
    end
  end
end

% Get rid of point regions
keepr = abs(R(:,1) - R(:,2)) > eps;
R = R(keepr,:);

num_regions = size(R,1);

% Create lookup table of which regions make up a time window
time2sreg_tab = cell(T,1);
for i = 1:T
  time2sreg_tab{i} = [];
end

for tidx = 1:T
  tj = t(tidx);
  Rjinds = [];
  for r = 1:num_regions
    if R(r,1) >= (tj - L) && R(r,2) <= (tj + L)
      Rjinds = [Rjinds, r]; %#ok
    end
  end
  time2sreg_tab{tidx} = Rjinds;
end

% Create lookup table mapping each region to the time indices that use this
% region
sreg2time_tab = cell(num_regions,1);
for i = 1:num_regions
  sreg2time_tab{i} = [];
end

for r = 1:num_regions
  tjs = [];
  for tidx = 1:T
    tj = t(tidx);
    if R(r,1) >= (tj - L) && R(r,2) <= (tj + L)
      tjs = [tjs, tidx]; %#ok
    end
  end
  sreg2time_tab{r} = tjs;
end
