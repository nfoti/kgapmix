function [num_subreg, measure, doc_subreg] = get_subreg(time_vec, triang_wt)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function takes as input a vector of times
% It also takes a functions that calculates the weight of a triangle whose tip
%   is a y-position y, the weight being calculated wrt some measure
%
% Consider a set of inverted right-angled isosceles trianges whose tips are at
% (t,0), where t belongs to time_vec. This function calculates the minimum 
% number of disjoint subregions whose union forms this set of triangles.
% It returns:
%  num_subreg: The number of disjoint subregions
%  measure   : The weights of each of these subregions
%  doc_subreg: The subregions which constitute each triangle
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

num_times = length(time_vec);

num_subreg = num_times * (num_times + 1)/2;

subreg_cnt = 0;

for i = 1:num_times
  doc_subreg{i} = [];
end;

for i = 1:num_times
  for j= 1:i
    for k = j:(j-1+num_times-(i-1))
      doc_subreg{k} = [doc_subreg{k},subreg_cnt+j];
    end;
    x(subreg_cnt+j)  = (-time_vec(j) + time_vec(((j-1)+num_times-(i-1))))/2;

    measure(subreg_cnt+j) = triang_wt(x(subreg_cnt+j));

    if i ~= 1
      for k = [max(1,j-1):min(i-1,j)]
        measure(subreg_cnt+j) = measure(subreg_cnt+j) - ...
                              triang_wt(x(k+subreg_cnt-(i-1)));
      end;

      if(j>1 & j < i)
        measure(subreg_cnt+j) = measure(subreg_cnt+j) + ...
                              triang_wt(x(j-1+subreg_cnt-(2*i-3)));
      end;
    end;
  end;
  subreg_cnt = subreg_cnt + i;
end;
