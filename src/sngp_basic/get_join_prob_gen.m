% VALIDATE BEFORE USING!!!
function p = get_join_prob(var_x, var_theta, measure, clstr_subrg, v_s, X, c,...
                           num_pts, sum_pts, indx)

% We are given the set of data values and cluster assignments. We want the
% probability that the points indexed by index joined the restaurant the way
% they did
% NOTE: all indexed points must belong to the same document !!!

  p = 0;
  clstrs = c(indx);
  x_list = X(indx);
  cl_id  = unique(clstrs);
  
  cl_num = length(cl_id);
  cl_cnt = zeros(cl_num,1);
  for i = 1:cl_num
    num_pts(cl_id(i)) = num_pts(cl_id(i)) - sum(clstrs == cl_id(i));
    sum_pts(cl_id(i)) = sum_pts(cl_id(i)) - sum(X(clstrs == cl_id(i)));
  end;

  for i=[1:length(indx)]
    c_curr = clstrs(i);

    mean_theta_post = sum_pts(c_curr) ./ (num_pts(c_curr) + var_x/var_theta);
    var_theta_post  = var_x*var_theta ./ (var_x + num_pts(c_curr).*var_theta);

    mean_post =(x_list(i).*var_theta_post + mean_theta_post.*var_x)./ ...
                (var_theta_post + var_x);
    var_post  = var_theta_post .* var_x ./ (var_theta_post + var_x);

    p = p + 0.5 * (log(var_post) -log(var_theta_post) + ...
                 mean_theta_post^2/var_theta_post - mean_post^2/var_post);

    if(num_pts(c_curr) == 0)
      num = measure(clstr_subrg(c_curr))
    else
      num = num_pts(c_curr);
    end;
    p = p + log(num) + log(v_s(clstr_subrg(c_curr)));
    num_pts(c_curr) = num_pts(c_curr) + 1;
    sum_pts(c_curr) = sum_pts(c_curr) + x_list(i);

  end;
