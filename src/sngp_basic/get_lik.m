function lik = get_lik(x_curr, num_pts, sum_pts, var_theta, var_x)

% Calculates the likelihood that a point joins a cluster with given suff.
% statistics

  mean_theta_post = sum_pts ./ (num_pts + var_x/var_theta);
  var_theta_post  = var_x*var_theta ./ (var_x + num_pts*var_theta);
  lik = log(normpdf(x_curr, mean_theta_post, sqrt(var_theta_post+var_x)));
